# Params
layer = 'RNA_counts'

out_files = []

for genes_name, genes in genes_dict.items():
    for reference in (references := scrnaseq_processed_data_dir.iterdir()):
        if 'lung' not in reference.stem or 'matched' not in reference.stem:
            continue

        reference_name = reference.stem
        reference_dir = seurat_to_h5_dir / reference_name
        reference_is_done = reference_dir / '.done'

        out_file = results_dir / f'embed_panel_restricted_genes_and_samples_scrnaseq/{genes_name}/{reference_name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
        out_files.append(out_file)

        rule:
            name: f'embed_panel_restricted_genes_and_samples_scrnaseq/{genes_name}/{reference_name}'
            input:
                reference_is_done=reference_is_done
            output:
                out_file=out_file,
            params:
                reference=reference_dir,
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
                mem='100GB',
                runtime='8h',
                # slurm_partition = "gpu",
                # slurm_extra = '--gres=gpu:1',
            conda:
                "spatial"
            shell:
                """
                mkdir -p "$(dirname {output.out_file})"

                python -u workflow/scripts/scRNAseq/embed_panel_scrnaseq.py \
                    --reference {params.reference} \
                    --out_file {output.out_file} \
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


rule embed_panel_restricted_genes_and_samples_scrnaseq_all:
    input:
        out_files