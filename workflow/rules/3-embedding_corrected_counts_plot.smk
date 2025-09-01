params_product = list(product(normalisations, layers, references, methods, levels))

out_files_panel = []

for correction_method in correction_methods:
    if correction_method =='raw':
        continue
    for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
        if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
            continue
        for condition in (conditions := segmentation.iterdir()): 
            for panel in (panels := condition.iterdir()):
                for normalisation, layer, reference, method, level in params_product:
                    if level not in CONDITIONS_LEVELS[condition.stem]:
                        continue
                    if reference not in CONDITIONS_REFERENCES[condition.stem]:
                        continue

                    # input embedding file (doesn't depend on ref,method or color loops but more readable to have here)
                    k = (segmentation.stem,condition.stem,panel.stem,normalisation)
                    name = '/'.join(k)                                                      
                    embed_file = results_dir / f'{correction_method}_embed_panel/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet'

                    # no need to plot panel for panel level UMAPs
                    if level == 'panel':
                        continue
                    
                    # no need to plot sample coloring for every param combination
                    if level == 'sample' and (reference != references[0] or method != methods[0]):
                        continue
                                                                                                                                                    # _{layer}
                    out_file = figures_dir / f"{correction_method}_embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}_{reference}_{method}_{level}.{extension}"
                    out_files_panel.append(out_file)

                    rule:
                        name: f'{correction_method}_embed_panel_plot/{name}/umap_{reference}_{method}_{level}_{layer}'
                        input:
                            panel=panel,
                            embed_file=embed_file,
                        output:
                            out_file=out_file,
                        params:
                            cell_type_annotation_dir=cell_type_annotation_dir,
                            normalisation=normalisation,
                            reference=reference,
                            method=method,
                            color=level,
                            cell_type_palette=cell_type_palette,
                            panel_palette=panel_palette,
                            sample_palette=sample_palette,
                            s=s,
                            alpha=alpha,
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

                            python workflow/scripts/xenium/embed_panel_plot.py \
                            --panel {input.panel} \
                            --embed_file {input.embed_file} \
                            --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                            --normalisation {params.normalisation} \
                            --reference {params.reference} \
                            --method {params.method} \
                            --color {params.color} \
                            --out_file {output.out_file} \
                            --cell_type_palette {params.cell_type_palette} \
                            --panel_palette {params.panel_palette} \
                            --sample_palette {params.sample_palette} \
                            --s {params.s} \
                            --alpha {params.alpha} \
                            --points_only {params.points_only} \

                            echo "DONE"
                            """


rule embed_panel_corrected_counts_plot_all:
    input:
        out_files_panel
