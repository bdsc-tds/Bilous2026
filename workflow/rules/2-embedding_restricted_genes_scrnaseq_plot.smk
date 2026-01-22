s=3
layer = 'RNA_counts'

out_files_panel = []

for genes_name, genes in genes_dict.items():
    for reference in (references := scrnaseq_processed_data_dir.iterdir()):
        reference_name = reference.stem
        reference_dir = seurat_to_h5_dir / reference_name
        
        for color in colors:
            if color == 'Level2.1' and 'external' in reference_name:
                continue
            if 'melanoma' in reference_name or 'PDAC' in reference_name or 'CRC' in reference_name:
                continue

            # input embedding file (doesn't depend on ref,method or color loops but more readable to have here)
            embed_file = results_dir / f'embed_panel_restricted_genes_scrnaseq/{genes_name}/{reference_name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 

            # no need to plot panel for panel color UMAPs
            if color == 'panel':
                continue
            
            # no need to plot sample coloring for every param combination
            if color == 'sample':
                continue

            out_file = figures_dir / f"embed_panel_restricted_genes_scrnaseq/{genes_name}/{reference_name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}_{color}.{extension}"
            out_files_panel.append(out_file)

            rule:
                name: f'embed_panel_restricted_genes_scrnaseq_plot/{genes_name}/{reference_name}/umap_{layer}_{color}'
                input:
                    embed_file=embed_file,
                output:
                    out_file=out_file,
                params:
                    reference=reference_dir,
                    color=color,
                    cell_type_palette=cell_type_palette,
                    panel_palette=panel_palette,
                    sample_palette=sample_palette,
                    s=s,
                    alpha=alpha,
                    dpi=dpi,
                    points_only='--points_only' if points_only else '',
                threads: 1
                resources:
                    mem='30GB',
                    runtime='10m',
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file})"

                    python workflow/scripts/scRNAseq/embed_panel_scrnaseq_plot.py \
                    --embed_file {input.embed_file} \
                    --reference {params.reference} \
                    --color {params.color} \
                    --out_file {output.out_file} \
                    --cell_type_palette {params.cell_type_palette} \
                    --panel_palette {params.panel_palette} \
                    --sample_palette {params.sample_palette} \
                    --s {params.s} \
                    --alpha {params.alpha} \
                    --dpi {params.dpi} \
                    {params.points_only} \

                    echo "DONE"
                    """


rule embed_panel_restricted_genes_scrnaseq_plot_all:
    input:
        out_files_panel
