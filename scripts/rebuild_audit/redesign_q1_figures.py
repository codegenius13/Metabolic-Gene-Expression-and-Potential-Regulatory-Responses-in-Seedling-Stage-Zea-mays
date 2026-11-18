from pathlib import Path
import shutil

script = Path("generate_maize_figures.py")
source = script.read_text(encoding="utf-8")

backup = Path("generate_maize_figures_before_q1_redesign_20261112.py.bak")

if backup.exists():
    raise FileExistsError(f"Backup already exists: {backup}")

shutil.copy2(script, backup)
print(f"Backup created: {backup}")

# ------------------------------------------------------------
# REPLACEMENT 1: Q1 ENRICHMENT DOT PLOT
# ------------------------------------------------------------

start = source.index("def q1_enrichment_dotplot():")
end = source.index(
    "# ============================================================\n# Q1 FIGURE 2: PROCESS DIRECTION",
    start
)

new_dotplot = r'''
def q1_enrichment_dotplot():
    rows = read_tsv(
        Q1_TABLES / "Q1_metabolic_process_enrichment.tsv"
    )

    rows = [
        r for r in rows
        if number(r.get("p.adjust")) > 0
        and number(r.get("FoldEnrichment")) > 0
    ]

    if not rows:
        print("Skipping Q1 enrichment dot plot: no valid rows.")
        return

    rows.sort(
        key=lambda r: (
            STRESSES.index(r["Stress"])
            if r["Stress"] in STRESSES else 99,
            number(r["p.adjust"])
        )
    )

    n = len(rows)
    y = np.arange(n)

    fig = plt.figure(figsize=(16, max(7.5, 0.82 * n + 2.0)))

    gs = fig.add_gridspec(
        1, 3,
        width_ratios=[1.35, 4.0, 7.0],
        wspace=0.035
    )

    id_ax = fig.add_subplot(gs[0, 0])
    label_ax = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[0, 2])

    # Keep the three axes aligned.
    for text_ax in (id_ax, label_ax):
        text_ax.set_ylim(-0.5, n - 0.5)
        text_ax.invert_yaxis()
        text_ax.set_xlim(0, 1)
        text_ax.axis("off")

    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()

    # Separate GO identifiers from process names.
    for i, row in enumerate(rows):
        id_ax.text(
            0.98, i, row["ID"],
            ha="right", va="center",
            fontsize=10, fontweight="bold",
            color="#333333"
        )

        label_ax.text(
            0.02, i, row["Description"],
            ha="left", va="center",
            fontsize=10.5, fontweight="bold",
            color="#222222"
        )

        label_ax.text(
            0.02, i + 0.25, f'Stress: {row["Stress"]}',
            ha="left", va="center",
            fontsize=8.5,
            color=stress_color(row["Stress"])
        )

    # Plot enrichment values, using point size for gene count.
    for stress in STRESSES:
        indices = [
            i for i, row in enumerate(rows)
            if row["Stress"] == stress
        ]

        if not indices:
            continue

        ax.scatter(
            [number(rows[i]["FoldEnrichment"]) for i in indices],
            [y[i] for i in indices],
            s=[
                max(65, number(rows[i]["Count"], 0) * 9)
                for i in indices
            ],
            color=stress_color(stress),
            edgecolors="#222222",
            linewidths=0.8,
            alpha=0.92,
            label=stress,
            zorder=3
        )

    ax.set_xlabel(
        "Fold enrichment",
        fontsize=12, fontweight="bold",
        labelpad=10
    )

    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=10)

    ax.grid(
        axis="x",
        linestyle="--",
        linewidth=0.8,
        alpha=0.35,
        zorder=0
    )

    ax.set_axisbelow(True)

    ax.set_title(
        "GO biological-process enrichment among metabolic-enzyme DEGs",
        fontsize=14, fontweight="bold", pad=18
    )

    ax.legend(
        title="Stress condition",
        frameon=False,
        ncol=4,
        fontsize=9,
        title_fontsize=9,
        loc="lower right"
    )

    fig.text(
        0.5, 0.025,
        "Point size represents the number of contributing genes.",
        ha="center",
        fontsize=9,
        color="#555555"
    )

    fig.subplots_adjust(
        left=0.035, right=0.98,
        top=0.88, bottom=0.12
    )

    save_figure(
        fig, Q1_FIGURES,
        "Q1_GO_enrichment_dotplot"
    )


'''

source = source[:start] + new_dotplot + source[end:]


# ------------------------------------------------------------
# REPLACEMENT 2: Q1 EXPRESSION DIRECTION
# ------------------------------------------------------------

start = source.index("def q1_process_direction():")
end = source.index(
    "# ============================================================\n# Q1 FIGURE 3: SHARED PROCESS ENRICHMENT HEATMAP",
    start
)

new_direction = r'''
def q1_process_direction():
    rows = read_tsv(
        Q1_TABLES / "Q1_metabolic_process_direction.tsv"
    )

    rows = [
        r for r in rows
        if not math.isnan(number(r.get("Median_log2FC")))
    ]

    if not rows:
        print("Skipping Q1 process direction plot: no valid rows.")
        return

    rows.sort(
        key=lambda r: (
            STRESSES.index(r["Stress"])
            if r["Stress"] in STRESSES else 99,
            r["Description"]
        )
    )

    n = len(rows)
    y = np.arange(n)

    fig = plt.figure(figsize=(19, max(7.5, 0.88 * n + 2.0)))

    gs = fig.add_gridspec(
        1, 4,
        width_ratios=[1.35, 4.0, 7.0, 3.2],
        wspace=0.035
    )

    id_ax = fig.add_subplot(gs[0, 0])
    label_ax = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[0, 2])
    class_ax = fig.add_subplot(gs[0, 3])

    # Align the information columns and the bar plot.
    for text_ax in (id_ax, label_ax, class_ax):
        text_ax.set_ylim(-0.5, n - 0.5)
        text_ax.invert_yaxis()
        text_ax.axis("off")

    id_ax.set_xlim(0, 1)
    label_ax.set_xlim(0, 1)
    class_ax.set_xlim(0, 1)

    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()

    values = [
        number(row["Median_log2FC"])
        for row in rows
    ]

    max_value = max(abs(v) for v in values)
    limit = max_value * 1.38 if max_value > 0 else 1

    ax.set_xlim(-limit, limit)

    # Separate GO IDs, process names and stress labels.
    for i, row in enumerate(rows):
        id_ax.text(
            0.98, i, row["GO_ID"],
            ha="right", va="center",
            fontsize=9.5, fontweight="bold",
            color="#333333"
        )

        label_ax.text(
            0.02, i, row["Description"],
            ha="left", va="center",
            fontsize=10, fontweight="bold",
            color="#222222"
        )

        label_ax.text(
            0.02, i + 0.25, f'Stress: {row["Stress"]}',
            ha="left", va="center",
            fontsize=8.5,
            color=stress_color(row["Stress"])
        )

    # Horizontal diverging bars.
    for i, row in enumerate(rows):
        value = values[i]

        ax.barh(
            i,
            value,
            color=stress_color(row["Stress"]),
            edgecolor="#222222",
            linewidth=0.65,
            alpha=0.93,
            height=0.62,
            zorder=3
        )

        # Place numerical values outside the bar endpoints.
        offset = limit * 0.025

        if value >= 0:
            ax.text(
                value + offset, i,
                f"+{value:.2f}",
                ha="left", va="center",
                fontsize=9, fontweight="bold",
                color="#222222"
            )
        else:
            ax.text(
                value - offset, i,
                f"{value:.2f}",
                ha="right", va="center",
                fontsize=9, fontweight="bold",
                color="#222222"
            )

    ax.axvline(
        0, color="#222222",
        linewidth=1.3, zorder=2
    )

    ax.set_xlabel(
        "Median log2 fold change",
        fontsize=11, fontweight="bold",
        labelpad=10
    )

    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=9)

    ax.grid(
        axis="x",
        linestyle="--",
        linewidth=0.8,
        alpha=0.35,
        zorder=0
    )

    ax.set_axisbelow(True)

    ax.set_title(
        "Expression direction of enriched biological processes",
        fontsize=14, fontweight="bold", pad=18
    )

    # Explicitly display the classification and gene counts.
    for i, row in enumerate(rows):
        direction = row["Direction"]

        if direction == "Predominantly upregulated":
            label = "Predominantly upregulated"
            color = "#287A4B"
        elif direction == "Predominantly downregulated":
            label = "Predominantly downregulated"
            color = "#B44444"
        else:
            label = "Mixed"
            color = "#B7791F"

        up = int(row["Up_Count"])
        down = int(row["Down_Count"])

        class_ax.text(
            0.02, i - 0.07, label,
            ha="left", va="center",
            fontsize=9.5, fontweight="bold",
            color=color
        )

        class_ax.text(
            0.02, i + 0.19,
            f"{up} up / {down} down",
            ha="left", va="center",
            fontsize=8.5,
            color="#555555"
        )

    # Stress legend.
    handles = [
        plt.Rectangle(
            (0, 0), 1, 1,
            color=COLORS[s],
            label=s
        )
        for s in STRESSES
        if any(row["Stress"] == s for row in rows)
    ]

    ax.legend(
        handles=handles,
        title="Stress condition",
        frameon=False,
        ncol=4,
        fontsize=8.5,
        title_fontsize=8.5,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01)
    )

    fig.text(
        0.5, 0.025,
        "Bars show median log2 fold change. Positive values extend right; "
        "negative values extend left. Direction labels use the 70% "
        "expression-direction classification rule.",
        ha="center",
        fontsize=8.5,
        color="#555555"
    )

    fig.subplots_adjust(
        left=0.035, right=0.985,
        top=0.84, bottom=0.12
    )

    save_figure(
        fig, Q1_FIGURES,
        "Q1_metabolic_process_expression_direction"
    )


'''

source = source[:start] + new_direction + source[end:]

script.write_text(source, encoding="utf-8")

print("Successfully replaced only the two Q1 figure functions.")
print("The shared-process heatmap function was not modified.")
print("The Q2 figure functions were not modified.")
