out_files_panel = []

for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for normalisation in normalisations: 
                for layer in layers:
                    k = (segmentation.stem,condition.stem,panel.stem,normalisation)
                    name = '/'.join(k)
                    rule_name = '/'.join(k+(layer,))

                    out_file = results_dir / f'embed_panel/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                    out_files_panel.append(out_file)

                    rule:
                        name: f'embed_panel/{rule_name}'
                        input:
                            panel=panel,
                        output:
                            out_file=out_file,
                        params:
                            xenium_processed_data_dir=xenium_processed_data_dir,
                            normalisation=normalisation,
                            layer=layer,
                            n_comps=n_comps,
                            n_neighbors=n_neighbors,
                            metric=metric,
                            min_dist=min_dist,
                            min_counts=min_counts,
                            min_features=min_features,
                            max_counts=max_counts,
                            max_features=max_features,
                            min_cells=min_cells,
                        threads: 1
                        resources:
                            mem='100GB' if panel.stem == '5k' else '50GB',
                            # runtime='30m' if panel.stem == '5k' else '20m',
                            runtime='8h' if panel.stem == '5k' else '3h',
                            # slurm_partition = "gpu",
                            # slurm_extra = '--gres=gpu:1',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file})"

                            python workflow/scripts/xenium/embed_panel.py \
                                --panel {input.panel} \
                                --out_file {output.out_file} \
                                --xenium_processed_data_dir {params.xenium_processed_data_dir} \
                                --normalisation {params.normalisation} \
                                --layer {params.layer} \
                                --n_comps {params.n_comps} \
                                --n_neighbors {params.n_neighbors} \
                                --metric {params.metric} \
                                --min_dist {params.min_dist} \
                                --min_counts {params.min_counts} \
                                --min_features {params.min_features} \
                                --max_counts {params.max_counts} \
                                --max_features {params.max_features} \
                                --min_cells {params.min_cells} \
                                
                            echo "DONE"
                            """


out_files_condition = []
for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
        continue
    for condition in (conditions := segmentation.iterdir()):
        n_panels = len(list(condition.iterdir()))
        if n_panels == 1: # no point to do condition plot if there's only one panel
            continue 
        for normalisation in normalisations: 
            for layer in layers:
                k = (segmentation.stem,condition.stem,normalisation)
                name = '/'.join(k)
                rule_name = '/'.join(k+(layer,))

                out_file = results_dir / f'embed_condition/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                out_files_condition.append(out_file)

                rule:
                    name: f'embed_condition/{rule_name}'
                    input:
                        condition=condition,
                    output:
                        out_file=out_file,
                    params:
                        xenium_processed_data_dir=xenium_processed_data_dir,
                        normalisation=normalisation,
                        layer=layer,
                        n_comps=n_comps,
                        n_neighbors=n_neighbors,
                        metric=metric,
                        min_dist=min_dist,
                        min_counts=min_counts,
                        min_features=min_features,
                        max_counts=max_counts,
                        max_features=max_features,
                        min_cells=min_cells,
                    threads: 1
                    resources:
                        mem='100GB' if panel.stem == '5k' else '50GB',
                        runtime='30m' if panel.stem == '5k' else '20m',
                        slurm_partition = "gpu",
                        slurm_extra = '--gres=gpu:1'
                    conda:
                        "spatial"
                    shell:
                        """
                        mkdir -p "$(dirname {output.out_file})"

                        python workflow/scripts/xenium/embed_condition.py \
                            --condition {input.condition} \
                            --out_file {output.out_file} \
                            --xenium_processed_data_dir {params.xenium_processed_data_dir} \
                            --normalisation {params.normalisation} \
                            --layer {params.layer} \
                            --n_comps {params.n_comps} \
                            --n_neighbors {params.n_neighbors} \
                            --metric {params.metric} \
                            --min_dist {params.min_dist} \
                            --min_counts {params.min_counts} \
                            --min_features {params.min_features} \
                            --max_counts {params.max_counts} \
                            --max_features {params.max_features} \
                            --min_cells {params.min_cells} \

                        echo "DONE"
                        """

rule embed_panel_all:
    input:
        out_files_panel

rule embed_condition_all:
    input:
        out_files_condition