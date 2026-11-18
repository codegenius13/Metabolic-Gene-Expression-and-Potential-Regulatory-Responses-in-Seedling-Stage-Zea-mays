#!/usr/bin/env python3
# Master maize RNA-seq rebuild. Run with: python3 master_maize_rebuild.py
# Requires Python 3.9+, Rscript, and the R package DESeq2. No pandas required.

import csv, gzip, hashlib, json, math, shutil, subprocess, sys, time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/mnt/c/Users/Boluwatife Ikuerowo/Desktop/Temp")
DE = ROOT / "05_DE_analysis"
ANN = ROOT / "06_pathway/annotation"
Q1 = ROOT / "07_results/01_answer"
Q2 = ROOT / "07_results/02_answer"
RDS = DE / "maize_DESeq2.rds"
GAF = ANN / "maize.B73.AGPv4.aggregate.gaf.gz"
OBO = ANN / "go-basic.obo"
ENZ = ANN / "enzyme_function_GO_terms.tsv"
TF = ANN / "regulatory_candidate_mapping/PlantTFDB_TF_mappings.tsv"
STAGE = ROOT / "07_results/_master_rebuild_staging"
STRESSES = ["Cold", "Drought", "Heat", "Salinity"]
ROOTS = {"GO:0008152","GO:0044238","GO:0044281","GO:0005975",
         "GO:0006082","GO:0006091","GO:0006520","GO:0006629"}
PADJ, LFC, MIN_GS, MAX_GS, DIRECTION = 0.05, 1.0, 5, 5000, 0.70

def say(s): print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)
def need(paths):
    missing = [str(p) for p in paths if not p.is_file()]
    if missing: raise FileNotFoundError("Missing required source files:\n" + "\n".join(missing))
def readtsv(p):
    with p.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))
def writetsv(p, rows, cols):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows: w.writerow({c: "" if r.get(c) is None else r.get(c) for c in cols})
def norm(x):
    x = (x or "").strip()
    return x[5:] if x.startswith("gene:") else x
def unique(xs, label):
    c = Counter(xs); dup = [x for x,n in c.items() if n > 1]
    if dup: raise ValueError(f"{label}: duplicate IDs after normalization, e.g. {dup[:8]}")
def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024), b""): h.update(b)
    return h.hexdigest()

def ontology():
    terms={}; cur=None
    with OBO.open(encoding="utf-8") as f:
        for line in f:
            s=line.rstrip()
            if s=="[Term]":
                if cur and cur.get("id"): terms[cur["id"]]=cur
                cur={"id":"","name":"","ns":"","alt":[],"parents":[],"obsolete":False,"replaced":[]}
            elif s.startswith("[") and s.endswith("]"):
                if cur and cur.get("id"): terms[cur["id"]]=cur
                cur=None
            elif cur is not None:
                if s.startswith("id: "): cur["id"]=s[4:]
                elif s.startswith("name: "): cur["name"]=s[6:]
                elif s.startswith("namespace: "): cur["ns"]=s[11:]
                elif s.startswith("alt_id: "): cur["alt"].append(s[8:])
                elif s.startswith("is_a: "): cur["parents"].append(s[6:].split()[0])
                elif s.startswith("relationship: part_of "): cur["parents"].append(s.split()[2])
                elif s=="is_obsolete: true": cur["obsolete"]=True
                elif s.startswith("replaced_by: "): cur["replaced"].append(s[13:].split()[0])
    if cur and cur.get("id"): terms[cur["id"]]=cur
    alt={}
    for gid,t in terms.items():
        for a in t["alt"]:
            if a in alt and alt[a]!=gid: raise ValueError("Conflicting GO alt_id: "+a)
            alt[a]=gid
    def resolve(g):
        g=alt.get(g,g); t=terms.get(g)
        if not t: return None
        if not t["obsolete"]: return g
        r=t["replaced"]
        if len(r)==1 and r[0] in terms and not terms[r[0]]["obsolete"]: return r[0]
        return None
    parents=defaultdict(set)
    for gid,t in terms.items():
        if not t["obsolete"]:
            for p in t["parents"]:
                rp=resolve(p)
                if rp: parents[gid].add(rp)
    memo={}
    def metabolic(g, trail=None):
        if g in ROOTS: return True
        if g in memo: return memo[g]
        trail=set() if trail is None else trail
        if g in trail: return False
        trail.add(g); v=any(metabolic(p,trail) for p in parents[g]); trail.remove(g)
        memo[g]=v; return v
    bp={g for g,t in terms.items() if not t["obsolete"] and t["ns"]=="biological_process" and metabolic(g)}
    names={g:t["name"] for g,t in terms.items() if not t["obsolete"]}
    return terms, resolve, bp, names, alt

def enzyme_terms(resolve):
    out=set()
    with ENZ.open(encoding="utf-8-sig") as f:
        for row in csv.reader(f, delimiter="\t"):
            for cell in row:
                for token in cell.replace(",", " ").replace(";", " ").split():
                    if token.startswith("GO:"):
                        r=resolve(token)
                        if r: out.add(r)
    if not out: raise ValueError("No GO IDs found in enzyme_function_GO_terms.tsv")
    return out

def annotations(resolve, bp_terms, enz):
    mf=defaultdict(set); bp=defaultdict(set); termgenes=defaultdict(set); audit=Counter()
    with gzip.open(GAF, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("!"): continue
            x=line.rstrip("\n").split("\t")
            if len(x)<9: audit["malformed"]+=1; continue
            gene,qual,go,aspect=x[1].strip(),x[3].strip(),x[4].strip(),x[8].strip()
            if not gene.startswith("Zm00001d"): continue
            audit["maize_records"]+=1
            if "NOT" in qual.split("|"): audit["NOT_excluded"]+=1; continue
            # Strict direct-evidence policy: accept only blank GAF qualifiers.
            if qual: audit["nonblank_qualifier_excluded"]+=1; continue
            go=resolve(go)
            if not go: audit["unresolved_or_obsolete_GO_excluded"]+=1; continue
            if aspect=="F" and go in enz: mf[gene].add(go); audit["enzyme_MF_accepted"]+=1
            if aspect=="P" and go in bp_terms: bp[gene].add(go); termgenes[go].add(gene); audit["metabolic_BP_accepted"]+=1
    eligible=set(mf)&set(bp)
    return mf,bp,termgenes,eligible,audit

def run_r(stage):
    rs=shutil.which("Rscript")
    if not rs: raise RuntimeError("Rscript unavailable in this Ubuntu environment.")
    rfile=stage/"extract_from_rds.R"
    rfile.write_text("""
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly=TRUE); dds <- readRDS(a[1]); out <- a[2]
if (!inherits(dds,"DESeqDataSet")) stop("RDS is not a DESeqDataSet")
if (paste(deparse(design(dds)), collapse = "") != "~condition") {
    stop("Unexpected design: ", paste(deparse(design(dds)), collapse = ""))
}
if (nrow(dds)!=25968 || ncol(dds)!=10) stop("Unexpected RDS dimensions")
if (anyDuplicated(rownames(dds)) || anyDuplicated(colnames(dds))) stop("Duplicate IDs")
if (!all(c("Control","Cold","Drought","Heat","Salinity") %in% as.character(colData(dds)$condition))) stop("Missing condition")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.table(data.frame(sample=rownames(colData(dds)),condition=as.character(colData(dds)$condition)),
 file=file.path(out,"sample_metadata.tsv"),sep="\\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
for (s in c("Cold","Drought","Heat","Salinity")) {
 z<-as.data.frame(results(dds,contrast=c("condition",s,"Control"),alpha=0.05)); z$Geneid<-rownames(z)
 z<-z[,c("Geneid",setdiff(names(z),"Geneid")),drop=FALSE]; z<-z[order(z$padj,na.last=TRUE),]
 write.table(z,file=file.path(out,paste0(s,"_vs_Control_all_results.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
}
""", encoding="utf-8")
    subprocess.run([rs,str(rfile),str(RDS),str(stage)],check=True)

def bh(ps):
    n=len(ps); order=sorted(range(n),key=lambda i:ps[i]); out=[1.0]*n; run=1.0
    for j in range(n-1,-1,-1):
        i=order[j]; run=min(run,ps[i]*n/(j+1)); out[i]=min(1.0,run)
    return out
def lc(n,k):
    if k<0 or k>n: return float("-inf")
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
def hyper_sf(km1,N,K,n):
    lo=max(km1+1,0,n-(N-K)); hi=min(K,n)
    if lo>hi: return 0.0
    den=lc(N,n)
    return min(1.0,sum(math.exp(lc(K,x)+lc(N-K,n-x)-den) for x in range(lo,hi+1)))

def main():
    need([RDS,GAF,OBO,ENZ,TF])
    if not (ROOT/"07_results").is_dir(): raise FileNotFoundError("07_results folder missing")
    if STAGE.exists(): shutil.rmtree(STAGE)
    for p in [STAGE/"de",STAGE/"ann",STAGE/"q1/data",STAGE/"q1/tables",
              STAGE/"q1/figures",STAGE/"q1/logs",STAGE/"q2/data",
              STAGE/"q2/tables",STAGE/"q2/figures",STAGE/"q2/logs"]:
        p.mkdir(parents=True,exist_ok=True)

    say("1/7 Extracting all four contrasts from the fitted DESeq2 RDS.")
    run_r(STAGE/"de")
    de={}
    for stress in STRESSES:
        rows=readtsv(STAGE/"de"/f"{stress}_vs_Control_all_results.tsv")
        ids=[norm(r["Geneid"]) for r in rows]; unique(ids,stress+" DESeq2 IDs")
        de[stress]={norm(r["Geneid"]):r for r in rows}
    universes=[set(de[s]) for s in STRESSES]
    if any(u!=universes[0] for u in universes[1:]): raise ValueError("Contrast gene universes differ")
    tested=universes[0]
    if len(tested)!=25968: raise ValueError("Expected 25,968 tested genes")

    say("2/7 Rebuilding GAF/OBO metabolic-enzyme and GO BP annotations.")
    terms,resolve,bp_terms,names,alt=ontology(); enz=enzyme_terms(resolve)
    mf,bp,termgenes,eligible,gaf_audit=annotations(resolve,bp_terms,enz)
    eligible_tested=eligible&tested
    if not eligible_tested: raise ValueError("No eligible genes overlap DESeq2 universe")
    annrows=[{"Gene_ID":g,"Enzyme_GO_Count":len(mf[g]),"Enzyme_GO_Terms":";".join(sorted(mf[g])),
              "Metabolic_BP_Count":len(bp[g]),"Metabolic_BP_Terms":";".join(sorted(bp[g]))} for g in sorted(eligible)]
    writetsv(STAGE/"ann/RefGen4_metabolic_enzyme_genes.tsv",annrows,
             ["Gene_ID","Enzyme_GO_Count","Enzyme_GO_Terms","Metabolic_BP_Count","Metabolic_BP_Terms"])
    writetsv(STAGE/"ann/RefGen4_metabolic_enzyme_genes.txt",[{"Gene_ID":g} for g in sorted(eligible)],["Gene_ID"])
    t2g={go:genes&eligible_tested for go,genes in termgenes.items()}
    writetsv(STAGE/"q1/data/metabolic_BP_TERM2GENE.tsv",
             [{"GO_ID":go,"Gene_ID":g} for go in sorted(termgenes) for g in sorted(termgenes[go])],
             ["GO_ID","Gene_ID"])
    writetsv(STAGE/"q1/data/metabolic_BP_TERM2NAME.tsv",
             [{"GO_ID":go,"Description":names.get(go,"")} for go in sorted(termgenes)],
             ["GO_ID","Description"])

    say("3/7 Regenerating metabolic-enzyme DEG sets.")
    candidates={}; candidate_rows={}; summary=[]
    for s in STRESSES:
        rr=[]
        for g,r in de[s].items():
            try: p=float(r["padj"]); l=float(r["log2FoldChange"])
            except (ValueError,TypeError): continue
            if g in eligible and math.isfinite(p) and math.isfinite(l) and p<=PADJ and abs(l)>=LFC:
                rr.append({"Geneid":g,"Gene_ID":g,**{k:r.get(k,"") for k in ["baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]}})
        rr.sort(key=lambda x:(float(x["padj"]),-abs(float(x["log2FoldChange"]))))
        candidates[s]={r["Gene_ID"] for r in rr}; candidate_rows[s]=rr
        cols=["Geneid","Gene_ID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]
        writetsv(STAGE/f"q1/data/{s}_metabolic_enzyme_DEGs.tsv",rr,cols)
        writetsv(STAGE/f"q2/data/{s}_metabolic_enzyme_candidates.tsv",rr,cols)
        summary.append({"Comparison":s+"_vs_Control","Total_tested_genes":len(de[s]),
                        "Metabolic_enzyme_DEGs":len(rr),
                        "Upregulated":sum(float(x["log2FoldChange"])>=1 for x in rr),
                        "Downregulated":sum(float(x["log2FoldChange"])<=-1 for x in rr)})
    writetsv(STAGE/"de/DEG_summary_metabolic_enzyme.tsv",summary,
             ["Comparison","Total_tested_genes","Metabolic_enzyme_DEGs","Upregulated","Downregulated"])
    say("4/7 Running per-stress GO enrichment and direction summaries.")

    bg = eligible_tested & (
        set().union(*t2g.values()) if t2g else set()
    )

    if not bg:
        raise ValueError("Empty Q1 enrichment background")

    enrichment = []
    directions = []

    for s in STRESSES:

        deg = candidates[s] & bg
        tests = []

        # Test EVERY eligible GO term, including terms with zero
        # DEG overlap, before applying BH correction.
        for go, genes in t2g.items():

            gs = genes & bg

            if len(gs) < MIN_GS or len(gs) > MAX_GS:
                continue

            overlap = gs & deg

            p = hyper_sf(
                len(overlap) - 1,
                len(bg),
                len(gs),
                len(deg)
            )

            tests.append({
                "ID": go,
                "Description": names.get(go, ""),
                "GeneRatio": f"{len(overlap)}/{len(deg)}",
                "BgRatio": f"{len(gs)}/{len(bg)}",
                "RichFactor": len(overlap) / len(gs),
                "FoldEnrichment": (
                    (len(overlap) / len(deg)) /
                    (len(gs) / len(bg))
                    if len(deg) > 0 else 0
                ),
                "pvalue": p,
                "geneID": "/".join(sorted(overlap)),
                "Count": len(overlap),
                "Stress": s
            })

        if not tests:
            continue

        # BH correction is performed over all tested GO terms
        # within this stress comparison.
        adjusted = bh([x["pvalue"] for x in tests])

        for x, adj in zip(tests, adjusted):

            x["p.adjust"] = adj
            x["qvalue"] = adj

            if adj <= PADJ and x["Count"] > 0:

                enrichment.append(x)

                vals = [
                    float(de[s][g]["log2FoldChange"])
                    for g in x["geneID"].split("/")
                ]

                up = sum(v >= LFC for v in vals)
                down = sum(v <= -LFC for v in vals)
                n = len(vals)

                direction = (
                    "Predominantly upregulated"
                    if up / n >= DIRECTION
                    else (
                        "Predominantly downregulated"
                        if down / n >= DIRECTION
                        else "Mixed"
                    )
                )

                import statistics

                directions.append({
                    "GO_ID": x["ID"],
                    "Description": x["Description"],
                    "Stress": s,
                    "DEG_Count": n,
                    "Up_Count": up,
                    "Down_Count": down,
                    "Up_Proportion": up / n,
                    "Down_Proportion": down / n,
                    "Median_log2FC": statistics.median(vals),
                    "Mean_log2FC": sum(vals) / n,
                    "Direction": direction
                })

    writetsv(
        STAGE / "q1/tables/Q1_metabolic_process_enrichment.tsv",
        enrichment,
        [
            "ID", "Description", "GeneRatio", "BgRatio",
            "RichFactor", "FoldEnrichment", "pvalue",
            "p.adjust", "qvalue", "geneID", "Count", "Stress"
        ]
    )

    writetsv(
        STAGE / "q1/tables/Q1_metabolic_process_direction.tsv",
        directions,
        [
            "GO_ID", "Description", "Stress", "DEG_Count",
            "Up_Count", "Down_Count", "Up_Proportion",
            "Down_Proportion", "Median_log2FC",
            "Mean_log2FC", "Direction"
        ]
    )

    writetsv(
        STAGE / "q1/tables/Q1_enrichment_background.tsv",
        [{"Gene_ID": g} for g in sorted(bg)],
        ["Gene_ID"]
    )

    bygo = defaultdict(list)

    for r in directions:
        bygo[r["GO_ID"]].append(r)

    overlaps = [
        {
            "GO_ID": g,
            "Description": v[0]["Description"],
            "Stress_Count": len(v),
            "Stresses": ";".join(x["Stress"] for x in v),
            "Directions": ";".join(
                x["Stress"] + ":" + x["Direction"]
                for x in v
            )
        }
        for g, v in bygo.items()
        if len(v) > 1
    ]

    writetsv(
        STAGE / "q1/tables/Q1_cross_stress_process_overlap.tsv",
        overlaps,
        [
            "GO_ID", "Description", "Stress_Count",
            "Stresses", "Directions"
        ]
    )

    say("5/7 Building Q2 candidates, GO links, and TF tables.")
    q2c=[]
    for s in STRESSES:
        for rank,r in enumerate(candidate_rows[s],1):
            q2c.append({"Stress":s,"Rank":rank,**r,"Enzyme_GO_Count":len(mf[r["Gene_ID"]]),
                        "Metabolic_BP_Count":len(bp[r["Gene_ID"]])})
    writetsv(STAGE/"q2/tables/Q2_metabolic_enzyme_candidates.tsv",q2c,
             ["Stress","Rank","Geneid","Gene_ID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","Enzyme_GO_Count","Metabolic_BP_Count"])
    q2go=[]
    for r in enrichment:
        for g in r["geneID"].split("/"):
            q2go.append({"Gene_ID":g,"Stress":r["Stress"],"GO_ID":r["ID"],
                         "GO_Description":r["Description"],"GO_adjusted_p_value":r["p.adjust"]})
    writetsv(STAGE/"q2/tables/Q2_candidate_to_GO_processes.tsv",q2go,
             ["Gene_ID","Stress","GO_ID","GO_Description","GO_adjusted_p_value"])
    oc=Counter(g for s in STRESSES for g in candidates[s])
    writetsv(STAGE/"q2/tables/Q2_candidate_cross_stress_overlap.tsv",
             [{"Stress_Count":n,"Unique_Genes":sum(v==n for v in oc.values())} for n in range(1,5)],
             ["Stress_Count","Unique_Genes"])

    tfrows=readtsv(TF)
    required={"TF_ID","Family","RefGen_v4_ID","Mapping_status"}
    if not tfrows or not required.issubset(tfrows[0]): raise ValueError("TF mapping schema is unexpected")
    fam=defaultdict(set); src=defaultdict(set); map_audit=Counter()
    for r in tfrows:
        status=(r["Mapping_status"] or "").strip().lower(); map_audit["records"]+=1; map_audit["status_"+status]+=1
        source=(r["TF_ID"] or "").strip(); target=norm(r["RefGen_v4_ID"]); family=(r["Family"] or "").strip()
        if status!="unique" or not source or not target or not family:
            map_audit["excluded_ambiguous_or_incomplete"]+=1; continue
        fam[target].add(family); src[target].add(source)
    tfgenes=set(fam)&tested
    tfsummary=[]; tfmatrix=[]; tffamily=[]; tfoverlap=[]
    for s in STRESSES:
        sig={}
        for g in tfgenes:
            r=de[s][g]
            try: p=float(r["padj"]); l=float(r["log2FoldChange"])
            except (ValueError,TypeError): continue
            if math.isfinite(p) and math.isfinite(l) and p<=PADJ and abs(l)>=LFC: sig[g]=(p,l)
        mapped=sum(g in de[s] for g in tfgenes)
        tfsummary.append({"Stress":s,"Mapped_TF_Genes_Tested":mapped,"Significant_TF_Genes":len(sig),
                          "Upregulated_Significant_TF_Genes":sum(l>=1 for p,l in sig.values()),
                          "Downregulated_Significant_TF_Genes":sum(l<=-1 for p,l in sig.values())})
        for g,(p,l) in sig.items():
            for f in sorted(fam[g]): tffamily.append({"Stress":s,"Gene_ID":g,"Family":f,"padj":p,"log2FoldChange":l})
            if g in candidates[s]:
                tfoverlap.append({"Stress":s,"Gene_ID":g,"TF_Families":";".join(sorted(fam[g])),"padj":p,"log2FoldChange":l})
    for g in sorted(tfgenes):
        row={"Gene_ID":g,"TF_Families":";".join(sorted(fam[g])),"PlantTFDB_Source_IDs":";".join(sorted(src[g]))}
        for s in STRESSES:
            r=de[s][g]; row[s+"_log2FoldChange"]=r["log2FoldChange"]; row[s+"_padj"]=r["padj"]
            try: row[s+"_significant"]="TRUE" if float(r["padj"])<=PADJ and abs(float(r["log2FoldChange"]))>=LFC else "FALSE"
            except (ValueError,TypeError): row[s+"_significant"]="FALSE"
        tfmatrix.append(row)
    writetsv(STAGE/"q2/tables/Q2_TF_stress_summary.tsv",tfsummary,
             ["Stress","Mapped_TF_Genes_Tested","Significant_TF_Genes","Upregulated_Significant_TF_Genes","Downregulated_Significant_TF_Genes"])
    writetsv(STAGE/"q2/tables/Q2_TF_family_summary.tsv",tffamily,["Stress","Gene_ID","Family","padj","log2FoldChange"])
    matrixcols=["Gene_ID"]+[s+"_"+f for s in STRESSES for f in ["log2FoldChange","padj","significant"]]+["TF_Families","PlantTFDB_Source_IDs"]
    writetsv(STAGE/"q2/tables/Q2_TF_cross_stress_expression_matrix.tsv",tfmatrix,matrixcols)
    writetsv(STAGE/"q2/tables/Q2_TF_metabolic_candidate_gene_overlap.tsv",tfoverlap,
             ["Stress","Gene_ID","TF_Families","padj","log2FoldChange"])

    say("6/7 Writing audits, manifest, and validating staged files.")
    writetsv(STAGE/"ann/GAF_annotation_audit.tsv",[{"Metric":k,"Count":v} for k,v in sorted(gaf_audit.items())],["Metric","Count"])
    writetsv(STAGE/"ann/GO_ontology_audit.tsv",
             [{"Metric":"GO_terms","Count":len(terms)},{"Metric":"GO_alt_IDs","Count":len(alt)},
              {"Metric":"Metabolic_BP_terms","Count":len(bp_terms)},{"Metric":"Enzyme_MF_terms","Count":len(enz)},
              {"Metric":"Eligible_metabolic_enzyme_genes","Count":len(eligible)},
              {"Metric":"Eligible_tested_genes","Count":len(eligible_tested)}],["Metric","Count"])
    writetsv(STAGE/"q2/logs/TF_mapping_audit.tsv",[{"Metric":k,"Count":v} for k,v in sorted(map_audit.items())],["Metric","Count"])
    # Hard checks: every Q2 GO contributor must be a candidate in that contrast.
    for r in q2go:
        if r["Gene_ID"] not in candidates[r["Stress"]]: raise ValueError("Invalid Q2 candidate-to-GO association")
    for r in enrichment:
        if not 0<=float(r["p.adjust"])<=1: raise ValueError("Invalid GO adjusted p-value")
    if any(len(de[s])!=25968 for s in STRESSES): raise ValueError("DESeq2 result row count mismatch")
    manifest={"created":time.strftime("%Y-%m-%d %H:%M:%S"),"root":str(ROOT),
              "thresholds":{"padj":PADJ,"abs_log2FC":LFC},
              "GO_policy":{"relationships":["is_a","part_of"],"qualifiers":"blank only; NOT and other qualifiers excluded",
                           "alt_ids":"resolved","obsolete":"single unambiguous replaced_by only"},
              "enrichment":{"background":"tested ∩ eligible metabolic-enzyme genes ∩ metabolic BP annotated genes",
                            "BH":"separately within each stress","min_gene_set":MIN_GS,"max_gene_set":MAX_GS},
              "counts":{"tested":len(tested),"eligible":len(eligible),"eligible_tested":len(eligible_tested),
                        "Q1_background":len(bg),"Q1_significant_term_stress":len(enrichment),
                        "Q2_candidate_GO_associations":len(q2go),"TF_target_genes_tested":len(tfgenes)},
              "source_sha256":{str(p.relative_to(ROOT)):digest(p) for p in [RDS,GAF,OBO,ENZ,TF]}}
    (STAGE/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    writetsv(STAGE/"validation_report.tsv",
             [{"Check":"DESeq2 contrasts","Status":"PASS","Detail":"Four contrasts; 25,968 genes each"},
              {"Check":"Gene ID normalization","Status":"PASS","Detail":"No collisions detected"},
              {"Check":"Metabolic annotation","Status":"PASS","Detail":str(len(eligible))+" eligible genes"},
              {"Check":"Q1 enrichment","Status":"PASS","Detail":str(len(enrichment))+" significant term-stress results"},
              {"Check":"Q2 candidate-GO links","Status":"PASS","Detail":str(len(q2go))+" associations"},
              {"Check":"TF mapping","Status":"PASS","Detail":str(len(tfgenes))+" unambiguous target genes tested"},
              {"Check":"Staging","Status":"PASS","Detail":"Active Q1/Q2 outputs untouched before validation"}],
             ["Check","Status","Detail"])

    say("7/7 Publishing validated outputs.")
    # Replace only generated folders/files; preserve scripts and root documents.
    tasks = [
        (
            STAGE / "de",
            DE,
            [
                "Cold_vs_Control_all_results.tsv",
                "Drought_vs_Control_all_results.tsv",
                "Heat_vs_Control_all_results.tsv",
                "Salinity_vs_Control_all_results.tsv",
                "sample_metadata.tsv",
                "R_sessionInfo.txt",
                "DEG_summary_metabolic_enzyme.tsv"
            ]
        ),
        (
            STAGE / "ann",
            ANN,
            [
                "RefGen4_metabolic_enzyme_genes.tsv",
                "RefGen4_metabolic_enzyme_genes.txt",
                "GAF_annotation_audit.tsv",
                "GO_ontology_audit.tsv"
            ]
        ),
        (
            STAGE / "q1",
            Q1,
            [
                "data",
                "tables",
                "logs"
            ]
        ),
        (
            STAGE / "q2",
            Q2,
            [
                "data",
                "tables",
                "logs"
            ]
        )
    ]
    rollback = ROOT / "07_results" / (
        "_master_rebuild_rollback_" +
        time.strftime("%Y%m%d_%H%M%S")
    )

    rollback.mkdir(parents=True, exist_ok=False)

    changed = []

    try:

        for srcroot, dstroot, names in tasks:

            for name in names:

                src = srcroot / name

                if not src.exists():
                    raise FileNotFoundError(
                        f"Expected staged output missing: {src}"
                    )

                dst = dstroot / name

                old = rollback / (
                    dstroot.name + "__" + name.replace("/", "_")
                )

                dstroot.mkdir(parents=True, exist_ok=True)

                # Preserve the original before replacing it.
                if dst.exists():

                    if dst.is_dir():
                        shutil.copytree(dst, old)
                    else:
                        shutil.copy2(dst, old)

                # Register the destination before attempting replacement,
                # so even a failed copy can be rolled back.
                changed.append((dst, old, dst.exists()))

                if dst.exists():
                    if dst.is_dir():
                        shutil.rmtree(dst)
                    else:
                        dst.unlink()

                if src.is_dir():
                    shutil.copytree(src, dst)
                else:
                    shutil.copy2(src, dst)

    except Exception:

        say("Publish error; attempting rollback.")

        for dst, old, existed_before in reversed(changed):

            if dst.exists():
                if dst.is_dir():
                    shutil.rmtree(dst)
                else:
                    dst.unlink()

            if existed_before and old.exists():

                if old.is_dir():
                    shutil.copytree(old, dst)
                else:
                    shutil.copy2(old, dst)

        say("Rollback attempt completed.")
        say(f"Rollback copies retained at: {rollback}")

        raise

    # Keep rollback copies until the new results have been reviewed.
    say(f"Previous outputs preserved at: {rollback}")

    shutil.copy2(
        STAGE / "manifest.json",
        ROOT / "07_results/rebuild_manifest.json"
    )

    shutil.copy2(
        STAGE / "validation_report.tsv",
        ROOT / "07_results/rebuild_validation_report.tsv"
    )

    shutil.rmtree(STAGE)

    say("SUCCESS. New validated outputs have replaced the active generated outputs.")
    say(f"Eligible metabolic-enzyme genes: {len(eligible)}")
    say(f"Q1 background: {len(bg)} | Q1 enriched term-stress results: {len(enrichment)}")
    say(f"Q2 candidate-GO associations: {len(q2go)} | TF target genes tested: {len(tfgenes)}")
    say("Review 07_results/rebuild_manifest.json and rebuild_validation_report.tsv.")
    say(f"Previous outputs are preserved at: {rollback}")

if __name__=="__main__":
    try: main()
    except Exception as e:
        say("STOPPED: "+str(e))
        say("Active Q1/Q2 outputs were not intentionally replaced before successful staging validation.")
        sys.exit(1)
