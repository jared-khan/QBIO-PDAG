# QBIO Public Data Analysis Group Final Project DESeq2 Code

This project aims to characterize differences in molecular signatures of breast cancer between young and old patients.

We used genomic and transcriptomic data from The Cancer Genome Atlas Program (TCGA), a large publicly available dataset provided by the National Cancer Institute and the National Human Genome Research Institute with 20,000 cancer and normal samples from 33 different cancer types. We also looked at proteomic data from the Clinical Proteomic Tumor Analysis Consortium (CPTAC), a smaller dataset also from the National Cancer Institute with proteomic data from a select number of the cancer samples from TCGA. We used this data to analyze and differentiate gene and protein expression in women with breast cancer in different age categories.

The overarching packages were TCGABiolinks for R, and cptac for Python.

For genomic analysis, the package maftools within TCGABiolinks was used to create the oncoplots based on mutation data from Mutation Annotation Format (MAF) data files. These oncoplots displayed the top mutated genes in young and old patients. Oncoplots were created using R and can be generated with Oncoplots.R.

Next, for a transcriptomic analysis, the gene expression for the top mutated genes was compared between the two age groups using boxplots. Gene expression data was acquired from the High Throughput Sequencing (HTSeq) counts files in the SummarizedExperiment package of TCGABiolinks. The HTSeq counts files are generated via HTSeqCounts.R.

Further transcriptomic analysis was performed through DESeq2 analysis from the DESeq2 and biomarRt R packages, which determined the most upregulated genes in terms of expression for young and old patients. The output for this can be recreated by running DESeq2Analysis.R.

SpearmanPlots.IPYNB contains the code to generate Spearman heatmaps, which were used to quantify the relationship
between the transcriptomics and proteomics for those genes. Additionall packages such as Numpy, Pandas, and Matplotlib were used to process and visualize the data.

Lastly, Kaplan-Meier curves using the survival package in TCGABiolinks were also generated for a few genes of interest from the oncoplots and the DESeq2 analysis. The KM plots specifically compared survival between young and old patients with high levels of expression for the specified genes. Patients with high expression were defined to have greater than the 75th percentile of the entire sample population, and this was determined through quartile analysis. The code for the KM plots are found in KaplainMeierCurves.R.
