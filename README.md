# Chromatinsight
An application to extract meaningful data from ChIP-seq experiments.

The code in Python represents the core part of the application, which uses the random forest model in the sklearn library (scikit-learn) to detect differential features between two groups of ChIP-seq experiments binarised using ChromHMM.

This software has been tested to detect the human sex dimorphism between two groups (males and females). The output can be further analysed using the R library associated, chromatinsight.tools https://github.com/ricolab/chromatinsight.tools . \

Input:
* The output text files of ChromHMM (one for each chromosome). At least ten samples in each group are needed to provide enough statistical power.
* A bed file with the genomic regions (such as TADs).

Output:
* A text file with the degree of dimorphism for each genomic region and each trial of the algorithm.
* A second text file (optional) with the degree of dimorphism of each genomic region after randomising the sample labels (important to display the random behaviour and calculating the FDR).

Requirements:\
Python 2.7\
pandas installed, see https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html \
sklearn installed, see https://scikit-learn.org/stable/install.html \
Quick way to install both:\
pip install pandas\
pip install sklearn
