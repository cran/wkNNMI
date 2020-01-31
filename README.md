# wkNNMI

wkNNMI is an R tool for the imputation of static and dynamic mixed-type data. A typical example of this kind of data are clinical registers containing subsequent screening visits for several patients.

This package implements an adaptive weighted k-nearest neighbours (wk-NN) imputation algorithm for clinical register data developed to explicitly handle missing values of continuous/ordinal/categorical and static/dynamic features conjointly. For each subject with missing data to be imputed, the method creates a feature vector constituted by the information collected over his/her first *window_size* time units of visits. This vector is used as sample in a k-nearest neighbours procedure, in order to select, among the other patients, the ones with the most similar temporal evolution of the disease over time. An *ad hoc* similarity metric was implemented for the sample comparison, capable of handling the different nature of the data, the presence of multiple missing values and include the cross-information among features.

## Installation
The package requires the following R packages to function correctly: *infotheo* and *foreach*.

#### Installation from CRAN

```r
install.packages("wkNNMI")
```

#### Installation from GitLab

wkNNMI is available on GitLab at https://gitlab.com/sysbiobig/wkNNMI

To install wkNNMI from GitLab, please use the following commands:
```r
library(devtools)
install_gitlab("sysbiobig/wkNNMI")
```

#### Installation from source package

The wkNNMI R source package can be downloaded at http://sysbiobig.dei.unipd.it/?q=wkNNMI

To install it from source, please use the following command:
```r
install.packages("wkNNMI_X.X.X.tar.gz", repos = NULL, type = "source")
```
where "X.X.X" indicates the package version.

## Getting started

The package contains two functions *impute.subject* and *impute.wknn*. Load the package and check the documentation of the two functions.

```r
library(wkNNMI)
?impute.subject
?impute.wknn
```

