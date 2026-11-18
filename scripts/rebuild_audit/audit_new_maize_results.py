#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict, Counter

ROOT = Path(__file__).resolve().parent
Q1 = ROOT / "07_results/01_answer/tables"
Q2 = ROOT / "07_results/02_answer/tables"

STRESSES = ["Cold", "Drought", "Heat", "Salinity"]
OUT = ROOT / "new_results_crosscheck.txt"


def read_table(path):
    if not path.exists():
        print(f"MISSING: {path}")
        return []

    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def integer(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def section(title):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


with OUT.open("w", encoding="utf-8") as report:
    import sys
    original_stdout = sys.stdout
    sys.stdout = report

    # --------------------------------------------------------
    # Q1: BACKGROUND
    # --------------------------------------------------------
    section("Q1 — BACKGROUND AND ANNOTATION")

    background = read_table(Q1 / "Q1_enrichment_background.tsv")
    print("Final enrichment background genes:", len({
        r.get("Gene_ID", "") for r in background if r.get("Gene_ID")
    }))

    # --------------------------------------------------------
    # Q1: ENRICHMENT
    # --------------------------------------------------------
    section("Q1 — ALL SIGNIFICANT ENRICHMENT RESULTS")

    enrichment = read_table(
        Q1 / "Q1_metabolic_process_enrichment.tsv"
    )

    print("Total term-stress rows:", len(enrichment))
    print("Distinct GO terms:", len({
        r.get("ID", "") for r in enrichment if r.get("ID")
    }))

    for stress in STRESSES:
        subset = [r for r in enrichment if r.get("Stress") == stress]
        print(f"\n{stress}: {len(subset)} enriched terms")

        for r in subset:
            print(
                f'  {r.get("ID")} | {r.get("Description")} | '
                f'GeneRatio={r.get("GeneRatio")} | '
                f'BgRatio={r.get("BgRatio")} | '
                f'RichFactor={r.get("RichFactor")} | '
                f'FoldEnrichment={r.get("FoldEnrichment")} | '
                f'p={r.get("pvalue")} | '
                f'padj={r.get("p.adjust")} | '
                f'q={r.get("qvalue")} | '
                f'Count={r.get("Count")} | '
                f'Genes={r.get("geneID")}'
            )

    # --------------------------------------------------------
    # Q1: PROCESS DIRECTION
    # --------------------------------------------------------
    section("Q1 — EXPRESSION DIRECTION FOR EACH ENRICHED PROCESS")

    directions = read_table(
        Q1 / "Q1_metabolic_process_direction.tsv"
    )

    for r in directions:
        print(
            f'{r.get("Stress")} | {r.get("GO_ID")} | '
            f'{r.get("Description")} | '
            f'DEG={r.get("DEG_Count")} | '
            f'Up={r.get("Up_Count")} | '
            f'Down={r.get("Down_Count")} | '
            f'Up%={r.get("Up_Proportion")} | '
            f'Down%={r.get("Down_Proportion")} | '
            f'Median={r.get("Median_log2FC")} | '
            f'Mean={r.get("Mean_log2FC")} | '
            f'Direction={r.get("Direction")}'
        )

    # --------------------------------------------------------
    # Q1: SHARED PROCESSES
    # --------------------------------------------------------
    section("Q1 — SHARED PROCESS SUMMARY")

    shared = read_table(
        Q1 / "Q1_cross_stress_process_overlap.tsv"
    )

    for r in shared:
        print(
            f'{r.get("GO_ID")} | {r.get("Description")} | '
            f'Stress count={r.get("Stress_Count")} | '
            f'Stresses={r.get("Stresses")} | '
            f'Directions={r.get("Directions")}'
        )

    # --------------------------------------------------------
    # Q2: METABOLIC CANDIDATES
    # --------------------------------------------------------
    section("Q2 — CANDIDATE COUNTS AND EXPRESSION DIRECTIONS")

    candidates = read_table(
        Q2 / "Q2_metabolic_enzyme_candidates.tsv"
    )

    all_candidate_genes = set()
    stress_gene_sets = {}

    for stress in STRESSES:
        subset = [r for r in candidates if r.get("Stress") == stress]
        genes = {r.get("Gene_ID") for r in subset if r.get("Gene_ID")}
        stress_gene_sets[stress] = genes
        all_candidate_genes.update(genes)

        up = sum(
            1 for r in subset
            if (number(r.get("log2FoldChange")) or 0) > 0
        )
        down = sum(
            1 for r in subset
            if (number(r.get("log2FoldChange")) or 0) < 0
        )

        print(
            f"{stress}: rows={len(subset)}, unique genes={len(genes)}, "
            f"up={up}, down={down}, total={up + down}"
        )

    print("\nTotal candidate-stress records:", len(candidates))
    print("Unique candidate genes across all stresses:", len(all_candidate_genes))

    frequency = Counter()
    for gene in all_candidate_genes:
        frequency[gene] = sum(
            gene in stress_gene_sets[s] for s in STRESSES
        )

    for n in range(1, 5):
        print(
            f"Genes occurring in exactly {n} stress condition(s): "
            f"{sum(v == n for v in frequency.values())}"
        )

    # --------------------------------------------------------
    # Q2: CANDIDATE CROSS-STRESS TABLE
    # --------------------------------------------------------
    section("Q2 — CANDIDATE CROSS-STRESS OVERLAP TABLE")

    overlap = read_table(
        Q2 / "Q2_candidate_cross_stress_overlap.tsv"
    )

    for r in overlap:
        print(
            f'{r.get("Stress_Count")} stresses: '
            f'{r.get("Unique_Genes")} genes'
        )

    # --------------------------------------------------------
    # Q2: TOP FIVE CANDIDATES
    # --------------------------------------------------------
    section("Q2 — TOP FIVE CANDIDATES PER STRESS")

    for stress in STRESSES:
        subset = [
            r for r in candidates
            if r.get("Stress") == stress
        ]
        subset.sort(key=lambda r: integer(r.get("Rank")))

        print(f"\n{stress}")
        for r in subset[:5]:
            print(
                f'Rank {r.get("Rank")} | {r.get("Gene_ID")} | '
                f'log2FC={r.get("log2FoldChange")} | '
                f'padj={r.get("padj")} | '
                f'baseMean={r.get("baseMean")} | '
                f'Enzyme GO count={r.get("Enzyme_GO_Count")} | '
                f'Metabolic BP count={r.get("Metabolic_BP_Count")}'
            )

    # --------------------------------------------------------
    # Q2: CANDIDATE-GO ASSOCIATIONS
    # --------------------------------------------------------
    section("Q2 — CANDIDATE-TO-GO ASSOCIATIONS")

    candidate_go = read_table(
        Q2 / "Q2_candidate_to_GO_processes.tsv"
    )

    print("Total candidate-GO association records:", len(candidate_go))

    for stress in STRESSES:
        subset = [r for r in candidate_go if r.get("Stress") == stress]
        print(
            f"{stress}: {len(subset)} associations; "
            f"{len({r.get('Gene_ID') for r in subset})} unique genes"
        )

    # Compare gene sets for GO terms occurring in multiple stresses.
    go_stress_genes = defaultdict(lambda: defaultdict(set))
    go_descriptions = {}

    for r in candidate_go:
        go = r.get("GO_ID", "")
        stress = r.get("Stress", "")
        gene = r.get("Gene_ID", "")
        if go and stress and gene:
            go_stress_genes[go][stress].add(gene)
            go_descriptions[go] = r.get("GO_Description", "")

    print("\nShared candidate genes for GO terms found in multiple stresses:")

    for go, by_stress in sorted(go_stress_genes.items()):
        present = [s for s in STRESSES if by_stress.get(s)]
        if len(present) < 2:
            continue

        union = set().union(*(by_stress[s] for s in present))
        shared_genes = set.intersection(*(by_stress[s] for s in present))

        print(
            f"\n{go} | {go_descriptions.get(go)} | "
            f"stresses={';'.join(present)}"
        )
        print("Unique genes across these stresses:", len(union))
        print("Genes present in every listed stress:", len(shared_genes))
        print("Shared genes:", ";".join(sorted(shared_genes)))

    # --------------------------------------------------------
    # Q2: TF SUMMARY
    # --------------------------------------------------------
    section("Q2 — SIGNIFICANT TF SUMMARY")

    tf_summary = read_table(
        Q2 / "Q2_TF_stress_summary.tsv"
    )

    for r in tf_summary:
        print(
            f'{r.get("Stress")} | tested={r.get("Mapped_TF_Genes_Tested")} | '
            f'significant={r.get("Significant_TF_Genes")} | '
            f'up={r.get("Upregulated_Significant_TF_Genes")} | '
            f'down={r.get("Downregulated_Significant_TF_Genes")}'
        )

    # --------------------------------------------------------
    # Q2: TF FAMILY COUNTS
    # --------------------------------------------------------
    section("Q2 — TF FAMILY COUNTS")

    tf_families = read_table(
        Q2 / "Q2_TF_family_summary.tsv"
    )

    family_genes_all = defaultdict(set)
    family_genes_stress = defaultdict(lambda: defaultdict(set))

    for r in tf_families:
        family = r.get("Family", "").strip()
        gene = r.get("Gene_ID", "").strip()
        stress = r.get("Stress", "").strip()

        if family and gene:
            family_genes_all[family].add(gene)
            if stress:
                family_genes_stress[stress][family].add(gene)

    print("Distinct TF families:", len(family_genes_all))
    print("Family-gene-stress records:", len(tf_families))

    print("\nOverall unique genes per family across all stress records:")
    for family, genes in sorted(
        family_genes_all.items(),
        key=lambda x: (-len(x[1]), x[0])
    ):
        print(f"{family}: {len(genes)}")

    for stress in STRESSES:
        print(f"\n{stress} — top TF families")
        ranked = sorted(
            family_genes_stress[stress].items(),
            key=lambda x: (-len(x[1]), x[0])
        )
        for family, genes in ranked[:15]:
            print(f"{family}: {len(genes)}")

    # --------------------------------------------------------
    # Q2: TF CROSS-STRESS MATRIX AND EXAMPLES
    # --------------------------------------------------------
    section("Q2 — TF MATRIX EXAMPLE GENES")

    tf_matrix = read_table(
        Q2 / "Q2_TF_cross_stress_expression_matrix.tsv"
    )

    examples = [
        "Zm00001d000112",
        "Zm00001d000045",
        "Zm00001d000303",
    ]

    for gene in examples:
        matches = [r for r in tf_matrix if r.get("Gene_ID") == gene]

        if not matches:
            print(f"{gene}: NOT FOUND")
            continue

        r = matches[0]
        print(f"\n{gene} | Families={r.get('TF_Families')}")

        for stress in STRESSES:
            print(
                f'  {stress}: log2FC={r.get(stress + "_log2FoldChange")} | '
                f'padj={r.get(stress + "_padj")} | '
                f'significant={r.get(stress + "_significant")}'
            )

    # --------------------------------------------------------
    # BASIC CONSISTENCY CHECKS
    # --------------------------------------------------------
    section("BASIC CONSISTENCY CHECKS")

    print("Q1 enrichment rows:", len(enrichment))
    print("Q1 direction rows:", len(directions))
    print("Q2 candidate-GO rows:", len(candidate_go))
    print("Q2 candidate rows:", len(candidates))
    print("Q2 TF matrix rows:", len(tf_matrix))
    print("Q2 TF family rows:", len(tf_families))

    print("\nAudit completed.")

    sys.stdout = original_stdout

print(f"Audit saved to: {OUT}")
