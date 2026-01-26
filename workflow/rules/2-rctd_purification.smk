from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# placeholders
xe_path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/xenium_paper/xenium/processed/segmentation/10x_15um/NSCLC/chuvio/0PSV/0PSV_2/std_seurat_objects/preprocessed_seurat.rds"
rctd_path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/xenium_paper/xenium/processed/segmentation/10x_15um/NSCLC/chuvio/0PSV/0PSV_2/cell_type_annotation/reference_based/matched_reference/rctd_class_aware/Level3/single_cell/output.rds"


rule purify_with_rctd:
    input:
        xe_path=xe_path,
        rctd_path=rctd_path
    output:
        post_processed_rctd_df=results_dir / f'puRCTD/post_processed_rctd_df.parquet',
        post_processed_rctd=results_dir / f'puRCTD/post_processed_rctd.rds', # feel free not to store
        spatial_neighborhood_scores=results_dir / f'puRCTD/spatial_neighborhood_scores.parquet', # can be used for cluster separation or certainty eval
        transcriptomics_neighborhood_scores=results_dir / f'puRCTD/transcriptomics_neighborhood_scores.parquet', # required for one of purification method (`score_based`) and be used for segmentation evaluation
        fully_purified_counts=results_dir / f'puRCTD/fully_purified_counts.parquet', 
        spot_class_purified_counts=results_dir / f'puRCTD/spot_class_counts.parquet',
        score_based_purified_counts=results_dir / f'puRCTD/score_based_purified_counts.parquet'
    script:
      "workflow/scripts/rctd_purification/post_rctd_scores_and_purification.R"
