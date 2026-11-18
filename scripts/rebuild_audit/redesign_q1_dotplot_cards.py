from pathlib import Path
import shutil

script = Path("generate_maize_figures.py")
source = script.read_text(encoding="utf-8")

# Create a backup without overwriting an existing backup.
backup = script.with_name("generate_maize_figures_before_bubble_cards.py.bak")
counter = 1

while backup.exists():
    backup = script.with_name(
        f"generate_maize_figures_before_bubble_cards_{counter}.py.bak"
    )
    counter += 1

shutil.copy2(script, backup)
print(f"Backup created: {backup}")

start = source.index("def q1_enrichment_dotplot():")
end = source.index(
    "# ============================================================\n# Q1 FIGURE 2: PROCESS DIRECTION",
    start
)

new_function = r'''
def q1_enrichment_dotplot():
    rows = read_tsv(
        Q1_TABLES / "Q1_metabolic_process_enrichment.tsv"
    )

    # Retain only significantly enriched processes.
    rows = [
        r for r in rows
        if number(r.get("p.adjust")) > 0
        and number(r.get("p.adjust")) <= 0.05
        and number(r.get("Count"), 0) > 0
    ]

    if not rows:
        print("Skipping Q1 enrichment bubble plot: no valid rows.")
        return

    # Order processes within each stress by adjusted p-value.
    for row in rows:
        row["_significance"] = -math.log10(
            max(number(row["p.adjust"]), 1e-300)
        )

    fig, axes = plt.subplots(
        2, 2,
        figsize=(16, 12),
        gridspec_kw={
            "hspace": 0.42,
            "wspace": 0.30
        }
    )

    axes = axes.flatten()

    # Consistent bubble-size mapping across all four panels.
    def bubble_size(count):
        return 65 + 11 * count

    # Same red colour family, independently normalized per stress.
    cmap = plt.get_cmap("Reds")

    for ax, stress in zip(axes, STRESSES):

        subset = [
            row for row in rows
            if row["Stress"] == stress
        ]

        ax.set_facecolor("#FFFFFF")

        # Make each stress panel look like a separate card.
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("#B8B8B8")
            spine.set_linewidth(1.1)

        ax.tick_params(
            axis="both",
            which="both",
            length=0,
            labelsize=9
        )

        ax.set_title(
            stress.upper() + " STRESS",
            fontsize=13,
            fontweight="bold",
            color="#8B1A1A",
            pad=14
        )

        if not subset:
            ax.text(
                0.5, 0.5,
                "No significantly enriched process",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=10,
                color="#777777"
            )
            ax.set_axis_off()
            continue

        subset.sort(
            key=lambda row: row["_significance"],
            reverse=True
        )

        significance = np.array([
            row["_significance"] for row in subset
        ])

        counts = np.array([
            number(row["Count"], 0) for row in subset
        ])

        y = np.arange(len(subset))

        # Each stress has its own significance range and normalization.
        local_max = max(float(np.max(significance)), 1.0)
        norm = plt.Normalize(vmin=0, vmax=local_max * 1.08)

        scatter = ax.scatter(
            significance,
            y,
            s=[bubble_size(count) for count in counts],
            c=significance,
            cmap=cmap,
            norm=norm,
            edgecolors="#651010",
            linewidths=0.9,
            alpha=0.98,
            zorder=3
        )

        # Biological-process names and GO identifiers remain outside
        # the bubbles, on separate lines.
        labels = [
            f'{row["Description"]}\n{row["ID"]}'
            for row in subset
        ]

        ax.set_yticks(y)
        ax.set_yticklabels(
            labels,
            fontsize=9,
            fontweight="medium"
        )
        ax.invert_yaxis()

        # Independent horizontal significance range for each stress.
        ax.set_xlim(0, local_max * 1.30)

        ax.set_xlabel(
            r"$-\log_{10}$(adjusted p-value)",
            fontsize=10,
            fontweight="bold",
            labelpad=8
        )

        ax.grid(
            axis="x",
            linestyle="--",
            linewidth=0.75,
            color="#D5D5D5",
            alpha=0.75,
            zorder=0
        )
        ax.set_axisbelow(True)

        # Separate red intensity gradient for each stress panel.
        colorbar = fig.colorbar(
            scatter,
            ax=ax,
            orientation="vertical",
            fraction=0.045,
            pad=0.025
        )

        colorbar.set_label(
            r"$-\log_{10}$(adjusted p-value)",
            fontsize=8.5,
            fontweight="bold"
        )
        colorbar.ax.tick_params(labelsize=8)
        colorbar.outline.set_linewidth(0.7)

        ax.margins(y=0.22)

    # Shared bubble-size legend: count is comparable across panels.
    legend_counts = [10, 30, 50]

    size_handles = [
        plt.scatter(
            [],
            [],
            s=bubble_size(count),
            color="#C73535",
            edgecolors="#651010",
            linewidths=0.8,
            alpha=0.95,
            label=str(count)
        )
        for count in legend_counts
    ]

    fig.legend(
        handles=size_handles,
        title="Number of contributing DEGs",
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=3,
        frameon=False,
        fontsize=9,
        title_fontsize=9
    )

    fig.suptitle(
        "Significant metabolic Gene Ontology biological-process "
        "enrichment across stress conditions",
        fontsize=16,
        fontweight="bold",
        y=0.985
    )

    fig.text(
        0.5, 0.005,
        "Bubble size represents the number of contributing DEGs. "
        "Red intensity represents enrichment significance. "
        "Each stress panel has its own significance scale. "
        "Only processes with adjusted p-value <= 0.05 are displayed.",
        ha="center",
        fontsize=8.5,
        color="#555555"
    )

    fig.subplots_adjust(
        left=0.29,
        right=0.94,
        top=0.91,
        bottom=0.12
    )

    save_figure(
        fig,
        Q1_FIGURES,
        "Q1_GO_enrichment_dotplot"
    )


'''

source = source[:start] + new_function + source[end:]
script.write_text(source, encoding="utf-8")

print("Updated q1_enrichment_dotplot() only.")
print("Direction figure, shared heatmap and Q2 functions were not changed.")
