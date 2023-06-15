# Cancer Prediction with AI

This project is my graduation thesis, which is a bioinformatics about colon cancer (CRC). The main goal of the project is to develop a colorectal cancer classification software which is cheaper, based on machine learning and classification methods. The other goal is finding new potential survival markers (biomarkers) with survival analysis. More information in the [report](docs/Report_CPWAI.pdf).

## How The Modules Work?

**1. Preprocessing:** This is an R project. It merges multiple NCBI GEO datasets into one single dataset. It also removes batch-effect to make values more indiscrete, and decreases the data amount with differential gene expression (DGE) analysis.
**2. Training & Testing (Classification):** This is a Python project. It uses the merged dataset to train models and test each one. It finds the best parameters for each classifier with Grid Search and gives different kinds of performance metrics as result.
**3. Survival Analysis:** This is an R project. It uses overall survival (OS) data to run survival fit. It finds genes which can be potential survival markers.

## How To Run Preprocessing (1_dataset_preprocesses)?

### 1_get_and_merge_datasets_vX.R

1. Set working directory as desired.
2. Set `gset_names` list with proper GEO series (dataset) names.
3. For each dataset, set the column name which includes MSS/MSI values in `mss_colnames` list. (Set in the same order with dataset names.)
4. Set `exprs_file_name` and `mdata_file_name` as desired. (Must be .csv files.)
5. Run all of the code. (RStudio is recommended, but optional.)

### 2_remove_batch_effect_vX.R

1. Set working directory as input directory.
2. Set `exprs_file_name` and `mdata_file_name` input file names, which will be the output files of `1_get_and_merge_datasets_vX.R`. (Must be .csv files.)
3. Set `corrected_exprs_file_name`, `pca_before_file_name` and `pca_after_file_name` output file names as desired. (`corrected_exprs_file_name` must be .csv, `pca_before_file_name` and `pca_after_file_name` must be .svg files.)
4. Run all of the code. (RStudio is recommended, but optional.)

### 3_differential_gene_expression_vX.R

1. Set working directory as input directory.
2. Set `exprs_file_name` and `mdata_file_name` input file names, which will be the output files of `2_remove_batch_effect_vX.R`. (Must be .csv files.)
3. Set `deg_exprs_file_name` and `up_and_down_table_file_name` output file names as desired. (Must be .csv files.)
4. Run all of the code. (RStudio is recommended, but optional.)

## How To Run Training & Testing (2_training_and_testing)?

1. Changing input configurations is not recommended. Necessary folders must be created if do not exist.
2. Copy merged metadata and expression datasets to input path.
3. Set shuffle split configurations as desired.
4. Set grid search configurations as desired.
5. Set the desired classifiers with their parameters. Set `param_grid` as `{}` for default parameters.
6. (Optional but recommended) Create virtual environment and activate it with (for Linux and MacOS):
    ```
    pip install virtualenv
    cd /path/to/2_training_and_testing
    python -m venv /path/to/new/virtual/environment
    source venv/bin/activate
    ```
7. Install dependencies from requirements.txt:
    ```
    pip install -r requirements.txt
    ```
8. Run the code:
    ```
    python main_vX.py
    ```

## How To Run Survival Analysis (3_survival_analysis)?

1. Set working directory as input directory.
2. Set `gset_name` GEO dataset name, which must include overall survival data.
3. Set `exprs_file_name`, `mdata_file_name` and `up_and_down_table_file_name` input file names. (Must be .csv files.)
4. Set `survival_mdata_file_name`, `survival_exprs_file_name`, `survival_up_and_down_table_file_name` and `survival_p_values_file_name` output file names as desired. (Must be .csv files.)
5. For `gset_name` dataset, set column names `os_time_column` and `os_event_column` which include overall survival time and event values.
6. Run all of the code. (RStudio is recommended, but optional.)
