This repository contains the analysis pipeline and figures presented in the paper downstream of the [xenium preprocessing pipeline](https://github.com/bdsc-tds/xenium_analysis_pipeline/tree/main), specifically the evaluation of `SPLIT`, `resolvi` and `ovrlpy` methods in terms of specificity, sensitivity, scib metrics, and cosine similarity with pseudobulk snRNAseq across platforms.

## Getting Started

### Prerequisites
*   **Snakemake:** For workflow management.
*   **Conda/Mamba:** For environment management (found in `workflow/envs`).

### Running the Pipeline
The workflow uses Snakemake profiles for execution on a SLURM cluster:

```bash
./run.slurm
```

Since the pipeline is very computationally intensive, by default all steps are commented out in `workflow/Snakefile`.
The include statement for a `.smk` file and corresponding rules in `rule all` can be uncommented one at a time to run the step of interest. 

## Repository Structure
*   **`data/`**: Directories for `xenium` and `scRNAseq` raw/processed datasets and associated metadata.
*   **`config/`**: Snakemake pipeline paths defined in `config.yml`.
*   **`workflow/`**: Snakemake logic, including `rules/` and `scripts/`. 
*   **`workflow/notebooks/`**: Drafts used to create the snakemake rules. Only few notebooks were not incorporated as snakemake rules:
    *   `0-palettes.ipynb`: Creates figure legends used throughout the paper. 
    *   `xenium/1-xenium_decontamination_comparison_logreg.ipynb`: Creates some subplots of figure 4. 
    *   `xenium/1-xenium_revision_figure_*.ipynb`: Creates some of the supplementary figures. 
*   **`figures/`**:  Figures outputs and associated .csv data.
*   **`figures_manuscript/`**: Figures symlinks into paper figures organized by figure number.
*   **`figures_manuscript_data/`**: Figures data symlinks into paper figures organized by figure number.
*   **`results/`**: Analysis outputs, including:
    *   `scib_metrics`: Integration benchmarking results.
    *   `contamination_metrics`: Specificity/precision calculations.
    *   `resolvi`: Results from the Resolvi signal enhancement method.


## Key Analysis Modules
*   **Contamination Metrics:** Evaluates spatial specificity and log-regression based precision/recall.
*   **Resolvi & Ovrlpy:** Scripts for evaluating signal integrity thresholds (0.5, 0.7) and corrected counts.
*   **Separability Metrics:** Analysis of how well different cell types are resolved in spatial vs. dissociated platforms.
