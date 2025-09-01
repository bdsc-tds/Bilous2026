params_product = list(product(normalisations, layers, references, methods, levels, [c for c in correction_methods if c != 'raw']))

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
                            for normalisation, layer, reference, method, level, correction_method in params_product:
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

                                if correction_method == "split_fully_purified":
                                    name_corrected = f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/split_fully_purified/'
                                    sample_corrected_counts_path = xenium_count_correction_dir / f"{name_corrected}/corrected_counts.h5"

                                else:
                                    if correction_method in ["resolvi", "resolvi_panel_use_batch=True", "resolvi_panel_use_batch=False"]:
                                        name_corrected = f'{name}/{mixture_k=}/{num_samples=}/'
                                    elif correction_method in ["resolvi_supervised", "resolvi_panel_supervised_use_batch=True", "resolvi_panel_supervised_use_batch=False"]:
                                        name_corrected = f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/{mixture_k=}/{num_samples=}'
                                    elif "ovrlpy" in correction_method:
                                        name_corrected = f'{name}'

                                    sample_corrected_counts_path = results_dir / f"{correction_method}/{name_corrected}/corrected_counts.h5"

                                output_dir = f'contamination_metrics_{name_params}/{genes_name}'

                                sample_normalised_counts = std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/{layer}.parquet'
                                sample_idx = std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/cells.parquet'
                                sample_annotation = cell_type_annotation_dir / f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet'
                                precomputed_ctj_markers = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_marker_genes.parquet'
                                precomputed_adata_obs = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_adata_obs.parquet'

                                out_file_df_ctj_marker_genes = results_dir /  f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_marker_genes.parquet'
                                out_file_df_diffexpr = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_diffexpr.parquet'
                                out_file_df_markers_rank_significance_diffexpr = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_markers_rank_significance_diffexpr.parquet'
                                out_file_summary_stats = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_summary_stats.json'
                                out_file_adata_obs = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_adata_obs.parquet'

                                if sample_corrected_counts_path.exists():

                                    out_files.extend([
                                        out_file_df_ctj_marker_genes,
                                        out_file_df_diffexpr,
                                        out_file_df_markers_rank_significance_diffexpr,
                                        out_file_summary_stats,
                                        out_file_adata_obs
                                        ])

                                    rule:
                                        name: f'{output_dir}_corrected_counts/{correction_method}/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                                        input:
                                            sample_corrected_counts_path=sample_corrected_counts_path,
                                            sample_dir=sample_dir,
                                            sample_normalised_counts=sample_normalised_counts,
                                            sample_idx=sample_idx,
                                            sample_annotation=sample_annotation,
                                            precomputed_ctj_markers=precomputed_ctj_markers,
                                            precomputed_adata_obs=precomputed_adata_obs,
                                        output:
                                            out_file_df_ctj_marker_genes=out_file_df_ctj_marker_genes,
                                            out_file_df_diffexpr=out_file_df_diffexpr,
                                            out_file_df_markers_rank_significance_diffexpr=out_file_df_markers_rank_significance_diffexpr,
                                            out_file_summary_stats=out_file_summary_stats,
                                            out_file_adata_obs=out_file_adata_obs,
                                        params:
                                            radius=radius,
                                            n_repeats=n_repeats,
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
                                                --sample_corrected_counts_path {input.sample_corrected_counts_path} \
                                                --sample_dir {input.sample_dir} \
                                                --sample_normalised_counts {input.sample_normalised_counts} \
                                                --sample_idx {input.sample_idx} \
                                                --sample_annotation {input.sample_annotation} \
                                                --precomputed_ctj_markers {input.precomputed_ctj_markers} \
                                                --precomputed_adata_obs {input.precomputed_adata_obs} \
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


rule contamination_metrics_diffexpr_corrected_counts_all:
    input:
        out_files
    output:
        touch(results_dir / f"{output_dir}_corrected_counts.done")


