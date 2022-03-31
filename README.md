# ABCstats

## Overview

ABCstats is an R package for metabolomics data transformation. Adaptive Box-Cox (ABC) transformation improves the data normality for statistical analysis.

## Installation

**Check if R package "devtools" is installed.**
```
# If not, install R package "devtools" first
install.packages("devtools")
```

**Install ABCstats**
```
devtools::install_github("Waddlessss/ABCstats")
```

## Usage

### Input

ABCstats takes feature intensity table (dataframe) as input, with features in row and samples in column by default. 

Specifically, the column names are sample names. The first row includes sample group names. The first column should be unique identifier for each metabolic feature. The input feature intensity table can be prepared in .csv format. 

An example input data table can also be viewed in R:
```
# View the example input data table
View(DemoData)
```

In the first row, you need to label samples by their group names.

After preparing the input feature intensity table, you are ready to read it to R.

```
# Read the input feature intensity table in .csv format.
inputTable = read.csv(filename)
# Using demo data as input
# inputTable = DemoData
```

### ABC transformation main function
```
# Sample normalization using ABC transformation main function
TransformedTable = ABCtransform(DemoData)
```

### Output
The function `ABCtransform` returns the transformed feature intensity table (data frame). An extra column named `lambda` will be added to specify the optimized lambda value.

More descriptions are provided in R.
```
# Find more descriptions of the ABC transformation output
?ABCtransform()
```

## Citation
To be updated.
