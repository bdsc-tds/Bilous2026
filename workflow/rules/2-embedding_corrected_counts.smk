raw_corrected_counts = True
params_product = list(product(normalisations, layers, references, methods, levels, [c for c in correction_methods if c != 'raw']))

out_files_panel = []

for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for normalisation, layer, reference, method, level, correction_method in params_product:
                if level not in CONDITIONS_LEVELS[condition.stem] or level =='sample':
                    continue
                if reference not in CONDITIONS_REFERENCES[condition.stem]:
                    continue

                k = (segmentation.stem,condition.stem,panel.stem)
                name = '/'.join(k)

                if correction_method == 'split_fully_purified':
                    panel_path = xenium_count_correction_dir / name
                else:
                    panel_path = results_dir / f'{correction_method}/{name}'     
                                                                            
                out_file = results_dir / f'{correction_method}_embed_panel/{name}/{normalisation}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                out_files_panel.append(out_file)

                rule:
                    name: f'{correction_method}_embed_panel/{name}/{normalisation}_{layer}_{reference}'
                    input:
                        count_correction_is_done=results_dir / f"{correction_method}.done"
                    output:
                        out_file=out_file,
                    params:
                        panel=panel_path,
                        normalisation=normalisation,
                        layer=layer,
                        reference=reference,
                        method=method,
                        level=level,
                        n_comps=n_comps,
                        n_neighbors=n_neighbors,
                        metric=metric,
                        min_dist=min_dist,
                        min_counts=min_counts,
                        min_features=min_features,
                        max_counts=max_counts,
                        max_features=max_features,
                        min_cells=min_cells,
                        num_samples=num_samples,
                        mixture_k=mixture_k,
                        xenium_count_correction_dir=xenium_count_correction_dir,
                        results_dir=results_dir,
                        correction_method=correction_method,
                        cell_type_annotation_dir = cell_type_annotation_dir,
                        annotation_normalisation = annotation_normalisation,
                        raw_corrected_counts='--raw_corrected_counts' if raw_corrected_counts else '',
                    threads: 1
                    resources:
                        mem='100GB' if panel.stem == '5k' else '50GB',
                        # runtime='30m' if panel.stem == '5k' else '20m',
                        runtime='8h',
                        # slurm_partition = "gpu",
                        # slurm_extra = '--gres=gpu:1',
                    conda:
                        "spatial"
                    shell:
                        """
                        mkdir -p "$(dirname {output.out_file})"

                        python workflow/scripts/xenium/embed_panel_corrected_counts.py \
                            --panel {params.panel} \
                            --out_file {output.out_file} \
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
                            --num_samples {params.num_samples} \
                            --mixture_k {params.mixture_k} \
                            --xenium_count_correction_dir {params.xenium_count_correction_dir} \
                            --results_dir {params.results_dir} \
                            --correction_method {params.correction_method} \
                            --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                            --annotation_normalisation {params.annotation_normalisation} \
                            --reference {params.reference} \
                            --method {params.method} \
                            --level {params.level} \
                            {params.raw_corrected_counts}
                            
                        echo "DONE"
                        """

rule embed_panel_corrected_counts_all:
    input:
        out_files_panel