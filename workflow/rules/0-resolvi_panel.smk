out_files_training = []
for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
                    
            k = (segmentation.stem,condition.stem,panel.stem,f'{mixture_k=}')
            name = '/'.join(k)

            out_dir_resolvi_model = results_dir / f'resolvi_panel_{use_batch=}/{name}/model/'
            out_file_resolvi_model = out_dir_resolvi_model / 'model.pt'
            out_files_training.append(out_file_resolvi_model)

            rule:
                name: f'resolvi_panel_training/{name}'
                input:
                    panel=panel,
                output:
                    out_file_resolvi_model = out_file_resolvi_model
                params:
                    xenium_processed_data_dir=xenium_processed_data_dir,
                    out_dir_resolvi_model=out_dir_resolvi_model,
                    min_counts=min_counts,
                    min_features=min_features,
                    max_counts=max_counts,
                    max_features=max_features,
                    min_cells=min_cells,
                    max_epochs=max_epochs,
                    mixture_k=mixture_k,
                    use_batch='--use_batch' if use_batch else '',
                threads: 1
                resources:
                    mem='80GB',# if panel.stem == '5k' else '10GB',
                    runtime='8h',
                    slurm_partition = "gpu",
                    slurm_extra = '--gres=gpu:1',
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {params.out_dir_resolvi_model})"

                    python workflow/scripts/xenium/resolvi_panel_training.py \
                    --panel {input.panel} \
                    --xenium_processed_data_dir {params.xenium_processed_data_dir} \
                    --out_dir_resolvi_model {params.out_dir_resolvi_model} \
                    --min_counts {params.min_counts} \
                    --min_features {params.min_features} \
                    --max_counts {params.max_counts} \
                    --max_features {params.max_features} \
                    --min_cells {params.min_cells} \
                    --mixture_k {params.mixture_k} \
                    --max_epochs {params.max_epochs} \
                    {params.use_batch} \
                    
                    echo "DONE"
                    """


out_files_inference = []
for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):

            k_model = (segmentation.stem,condition.stem,panel.stem,f'{mixture_k=}')
            k = k_model[:-1]
            name_model = '/'.join(k_model)
            name = '/'.join(k)

            dir_resolvi_model = results_dir / f'resolvi_panel_{use_batch=}/{name_model}/model/'
            out_dir = results_dir / f'resolvi_panel_{use_batch=}/{name}/'
            out_file = out_dir / 'inference.done'
            out_files_inference.append(out_file)

            rule:
                name: f'resolvi_panel_inference/{name}'
                input:
                    panel=panel,
                    file_resolvi_model = dir_resolvi_model / 'model.pt'
                output:
                    out_file_inference=touch(out_file),
                params:
                    xenium_processed_data_dir=xenium_processed_data_dir,
                    dir_resolvi_model=dir_resolvi_model,
                    results_dir=results_dir,
                    min_counts=min_counts,
                    min_features=min_features,
                    max_counts=max_counts,
                    max_features=max_features,
                    min_cells=min_cells,
                    num_samples=num_samples,
                    batch_size=batch_size,
                    mixture_k=mixture_k,
                    use_batch='--use_batch' if use_batch else '',
                threads: 1
                resources:
                    mem='500GB',# if panel.stem == '5k' else '10GB',
                    runtime='12h',
                    slurm_partition = "gpu",
                    slurm_extra = '--gres=gpu:1',
                conda:
                    "spatial"
                shell:
                    """
                    python workflow/scripts/xenium/resolvi_panel_inference.py \
                    --panel {input.panel} \
                    --xenium_processed_data_dir {params.xenium_processed_data_dir} \
                    --dir_resolvi_model {params.dir_resolvi_model} \
                    --results_dir {params.results_dir} \
                    --min_counts {params.min_counts} \
                    --min_features {params.min_features} \
                    --max_counts {params.max_counts} \
                    --max_features {params.max_features} \
                    --min_cells {params.min_cells} \
                    --mixture_k {params.mixture_k} \
                    --num_samples {params.num_samples} \
                    --batch_size {params.batch_size} \
                    {params.use_batch} \

                    echo "DONE"
                    """


rule resolvi_panel_training_all:
    input:
        out_files_training

rule resolvi_panel_inference_all:
    input:
        out_files_inference
    output:
        touch(results_dir / f"resolvi_{use_batch=}.done")