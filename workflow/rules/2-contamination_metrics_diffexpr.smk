params_product = list(product(normalisations, layers, references, methods, levels))

out_files = []

for genes_name, genes in genes_dict.items():
    for markers_mode in markers_modes:
        for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
            if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
                continue
            for condition in (conditions := segmentation.iterdir()): 
                for panel in (panels := condition.iterdir()):
                    if genes_name != 'all' and panel.name != '5k':
                        continue
                    for donor in (donors := panel.iterdir()):
                        for sample in (samples := donor.iterdir()):
                            for normalisation, layer, reference, method, level in params_product:
                                if level not in CONDITIONS_LEVELS[condition.stem]:
                                    continue
                                if reference not in CONDITIONS_REFERENCES[condition.stem]:
                                    continue

                                k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                                name = '/'.join(k)

                                name_params = f"{markers_mode}_{radius=}_{top_n=}"

                                if 'proseg' in segmentation.stem:
                                    k_proseg = ('proseg',condition.stem,panel.stem,donor.stem,sample.stem)
                                    name_proseg = '/'.join(k_proseg)
                                    sample_dir = xenium_processed_data_dir / f'{name_proseg}/raw_results'
                                else:
                                    sample_dir = xenium_processed_data_dir / f'{name}/normalised_results/outs'
                                    

                                output_dir = f'contamination_metrics_{name_params}/{genes_name}'

                                sample_normalised_counts = std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/{layer}.parquet'
                                sample_idx = std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/cells.parquet'
                                sample_annotation = cell_type_annotation_dir / f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet'

                                out_file_df_ctj_marker_genes = results_dir /  f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_marker_genes.parquet'
                                out_file_df_diffexpr = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_diffexpr.parquet'
                                out_file_df_markers_rank_significance_diffexpr = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_markers_rank_significance_diffexpr.parquet'
                                out_file_summary_stats = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_summary_stats.json'
                                out_file_adata_obs = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_adata_obs.parquet'


                                out_files.extend([
                                    out_file_df_ctj_marker_genes,
                                    out_file_df_diffexpr,
                                    out_file_df_markers_rank_significance_diffexpr,
                                    out_file_summary_stats,
                                    out_file_adata_obs,
                                    ])

                                rule:
                                    name: f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                                    input:
                                        sample_dir=sample_dir,
                                        sample_normalised_counts=sample_normalised_counts,
                                        sample_idx=sample_idx,
                                        sample_annotation=sample_annotation,
                                    output:
                                        out_file_df_ctj_marker_genes=out_file_df_ctj_marker_genes,
                                        out_file_df_diffexpr=out_file_df_diffexpr,
                                        out_file_df_markers_rank_significance_diffexpr=out_file_df_markers_rank_significance_diffexpr,
                                        out_file_summary_stats=out_file_summary_stats,
                                        out_file_adata_obs=out_file_adata_obs,
                                    params:
                                        radius=radius,
                                        top_n=top_n,
                                        markers=markers_mode,
                                        min_counts=min_counts,
                                        min_features=min_features,
                                        max_counts=max_counts,
                                        max_features=max_features,
                                        min_cells=min_cells,
                                        genes=genes,
                                    threads: 1
                                    resources:
                                        mem='100GB',
                                        runtime='3h',
                                    conda:
                                        "spatial"
                                    shell:
                                        """
                                        mkdir -p "$(dirname {output.out_file_df_diffexpr})"

                                        python workflow/scripts/xenium/contamination_metrics_diffexpr_sample.py \
                                            --sample_dir {input.sample_dir} \
                                            --sample_normalised_counts {input.sample_normalised_counts} \
                                            --sample_idx {input.sample_idx} \
                                            --sample_annotation {input.sample_annotation} \
                                            --out_file_df_ctj_marker_genes {output.out_file_df_ctj_marker_genes} \
                                            --out_file_df_diffexpr {output.out_file_df_diffexpr} \
                                            --out_file_df_markers_rank_significance_diffexpr {output.out_file_df_markers_rank_significance_diffexpr} \
                                            --out_file_summary_stats {output.out_file_summary_stats} \
                                            --out_file_adata_obs {output.out_file_adata_obs} \
                                            --radius {params.radius} \
                                            --top_n {params.top_n} \
                                            --markers {params.markers} \
                                            --min_counts {params.min_counts} \
                                            --min_features {params.min_features} \
                                            --max_counts {params.max_counts} \
                                            --max_features {params.max_features} \
                                            --min_cells {params.min_cells} \
                                            --genes {params.genes} \

                                        echo "DONE"
                                        """


rule contamination_metrics_diffexpr_all:
    input:
        out_files
    output:
        touch(results_dir / f"{output_dir}.done")


