# Packgage: protqc

## Description
  The package protqc output Quality Control(QC) results of proteomics data for Quartet Project. The QC pipeline starts from the expression profiles at peptide/protein levels, and enables to calculate 6 metrics. A Total score is the geometric mean of the linearly normalized values of these metrics.<br />
  1.	***Number of features***: We expect as many proteins (mapped to gene symbols) as possible for downstreaming analyses.
  2.	***Missing percentage (%)***: Too many missing values interfere with comparability. This metric is calculated globally.
  3.	***Coefficient of variantion (CV, %)***: A CV value is calculated to indicate the dispersion within replicates feature by feature.
  4.	***Absolute Correlation***: Pearson correlation reflects overall reproducibility within replicates. We calculate correlation coefficients between each two replicates within each biological sample (D5, D6, F7, M8), and take the median as the final value for absolute correlation.
  5.	***Signal-to-Noise Ratio (SNR)***: SNR is established to characterize the ability of a platform or lab or batch, which is able to distinguish intrinsic differences among distinct biological sample groups (“signal”) from variations in technical replicates of the same sample group ("noise").
  6.	***Relative Correlation with Reference Datasets (RC)***: RC is used for assessment of quantitative consistency with the reference dataset at relative levels. For shotgun proteomics, quantitation at peptide levels is theoretically more reliable. Therefore, the reference dataset is established by benchmarking the relative expression values (log2FCs), for each peptide sequence of each sample pair (D5/D6, F7/D6, M8/D6), in historical datasets at peptide levels. We calculate relatively qualified (satisfied with thresholds of p < 0.05) log2FCs of the queried data, for overlapped peptides with the reference dataset, as the input for the assessment of quantitative consistency. Then RC value is Pearson correlation coefficient between the test dataset and the reference dataset.



## Depends
  Environment: R (>= 3.5.0)<br />
  Packages: dplyr,data.table,ggplot2,ggthemes,edgeR,limma,reshape2,psych

## Installation
```
devtools::install_github("QiaochuChen/protqc")
```

## Usage
  protqc::qc_conclusion(exp_path, meta_path, output_dir, plot)
  > You can also use other QC functions (start with 'qc') to get performance for your data in each metric.

## Examples
```
exp_path <- './test/input/example_data_for_test1.csv'
meta_path <- './test/input/example_meta_for_test1.csv'
output_dir <- './test/output/test1'
protqc::qc_conclusion(exp_path, meta_path, output_dir, plot=FALSE)
```

## Built-in Data
1. ./data/reference_dataset.rda<br />
   This file contains benchmarked historical data with relative quantities at peptide levels. For shot-gun proteomics, quantitation at peptide levels is theoretically more reliable. We use relatively quantitative expression values, for each sequence of each sample pair (D5/D6, F7/D6, M8/D6), as the reference for assessment of quantitative accuracy. 
   > For details in buiding the reference dataset, please refer to ... 
  
2. ./data/historical_data_genesymbols.rda<br />
   This file contains historical datasets (at protein levels, mapped to genes), which are used for calculating historical performances.

3. ./data/historical_data_peptides.rda<br />
   This file contains historical datasets (at peptide levels, unmapped sequences), which are used for calculating historical performances in quantitative accuracy at relative levels.

4. ./data/historical_meta.rda<br />
   This file contains historical metadata.

5. ./data/historical_qc.rda<br />
   ./data/historical_qc_norm.rda<br />
   ./data/historical_qc_stat.rda<br />
   These files are used for ranking scores of the testing data set in historical values.
   > These files are all output of the function ***qc_history()***.

## Examples for the input data
1. ./test/input/example_data_for_test1.csv & ./test/input/example_data_for_test2.csv<br />
   This file is an example for profiled data at protein levels and peptide levels (for these functions: ***input_data()***, ***qc_info()***, ***qc_snr()***, ***qc_cor()***, ***qc_conclusion()***), containing quantitative levels in each sample (replicate). The first two columns ("Type" and "Feature") are required to explain the feature type, which is either "Gene Symbol", or "Peptide Sequence") and the feature name. The missing values of data are replaced by *NA* or *0*.
   > If you only provide data at protein levels (refer to the file *example_data_for_test2.csv*), then the metric RC will not be calculated.

2. ./test/input/example_meta_for_test1.csv & ./test/input/example_meta_for_test2.csv<br />
   This file is an example for proteomic metadata. The columns "library" and "sample" are required.
   > The column names of profiled data (except the first two columns) and the column 'library' of metadata must be in one-to-one correspondence.
