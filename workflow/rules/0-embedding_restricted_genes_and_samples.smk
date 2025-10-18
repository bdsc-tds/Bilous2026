out_files_panel = []

for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
        continue
    for condition in (conditions := segmentation.iterdir()):
        for panel in (panels := condition.iterdir()):
            if panel.stem not in ['lung','5k']:
                continue
            for normalisation in normalisations: 
                for layer in layers:
                    k = (segmentation.stem,condition.stem,panel.stem,normalisation)
                    name = '/'.join(k)
                    rule_name = '/'.join(k+(layer,))

                    out_file = results_dir / f'embed_panel_restricted_genes_and_samples/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                    out_files_panel.append(out_file)

                    rule:
                        name: f'embed_panel_restricted_genes_and_samples/{rule_name}'
                        input:
                            panel=panel,
                        output:
                            out_file=out_file,
                        params:
                            xenium_processed_data_dir = xenium_processed_data_dir,
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
                            genes=nsclc_shared_genes,
                            samples=nsclc_shared_samples,
                        threads: 1
                        resources:
                            mem='100GB' if panel.stem == '5k' else '50GB',
                            # runtime='30m' if panel.stem == '5k' else '20m',
                            runtime='12h',
                            # slurm_partition = "gpu",
                            # slurm_extra = '--gres=gpu:1',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file})"

                            python -u workflow/scripts/xenium/embed_panel.py \
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
                                --genes {params.genes} \
                                --samples {params.samples} \
                                
                            echo "DONE"
                            """


rule embed_panel_restricted_genes_and_samples_all:
    input:
        out_files_panel