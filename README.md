# ProtQC

## Description
The package ProtQC visualizes Quality Control(QC) results for Quartet Project.

## Installation

```
devtools::install_github("chinese-quartet/ProtQC")
```

## Usage

ProtQC::plot_pca(expr_dt_path,meta_dt_path,output_dir)

ProtQC::plot_corr(expr_dt_path,meta_dt_path,output_dir)

ProtQC::table_conclusion(expr_dt_path,meta_dt_path,output_dir)


## Examples
```
ProtQC::plot_pca('./data/test_data_log2_example.csv','./data/test_metadata_example.csv','~/Desktop/')
```
```
ProtQC::plot_corr('./data/test_data_log2_example.csv','./data/test_metadata_example.csv','~/Desktop/')
```
```
ProtQC::table_conclusion('./data/test_data_log2_example.csv','./data/test_metadata_example.csv','~/Desktop/')
```

## Built-in Dataset

1. example_ref_dt.rds

> This file is an example for reference data set, which is a benchmarked proteomics profiled data set for Quality Control.

2. ref_snrcorr.rds

> This file is used in ranking scores of the testing data set in the historical values.

## Examples for the input data

1. data/input/test_data_log2_example.csv

This file is an example for proteomics profiled data (***expr_dt***, one of the input parameters for these functions: ***plot_pca()***, ***plot_corr()***,***table_conclusion()***), containing gene symbols of each protein and its quantitated expression level in each sample (replicate), and the missing values are replaced by *NA*. **Remember that each value is a base-2 logarithm of its original value.**

2. data/input/test_metadata_example.csv

This file is an example for proteomics metadata (***meta_dt***, one of the input parameters for these functions: ***plot_pca()***, ***plot_corr()***,***table_conclusion()***). 

Here is a quick view: 

| name                           | sample |
| :----------------------------- | ------ |
| B1_DDA_JNU_Lumos_Firmiana_D5_1 | D5     |
| B1_DDA_JNU_Lumos_Firmiana_D5_2 | D5     |
| B1_DDA_JNU_Lumos_Firmiana_D5_3 | D5     |
| B1_DDA_JNU_Lumos_Firmiana_D6_1 | D6     |
| B1_DDA_JNU_Lumos_Firmiana_D6_2 | D6     |
| B1_DDA_JNU_Lumos_Firmiana_D6_3 | D6     |
| B1_DDA_JNU_Lumos_Firmiana_F7_1 | F7     |
| B1_DDA_JNU_Lumos_Firmiana_F7_2 | F7     |
| B1_DDA_JNU_Lumos_Firmiana_F7_3 | F7     |
| B1_DDA_JNU_Lumos_Firmiana_M8_1 | M8     |
| B1_DDA_JNU_Lumos_Firmiana_M8_2 | M8     |
| B1_DDA_JNU_Lumos_Firmiana_M8_3 | M8     |

> Note: The column names of expr_dt and meta_dt$name must be in one-to-one correspondence.





