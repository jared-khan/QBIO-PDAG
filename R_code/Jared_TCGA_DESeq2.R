library(TCGAbiolinks)
library(DESeq2)
library(stringr)
library(dplyr)
library(tidyverse)
library(readxl)
library(biomaRt)

# REFERENCE: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# If testing on the 6 patients, get the patient barcodes
# clinical <- read.csv("tcga_brca_six_example_clinical.csv")
# barcodes <- as.character( clinical_file$barcode )

# Read in and transpose counts data
data = read.csv("Jeff_TCGA_counts_raw.csv", header=FALSE)
data = t(data)

# Takes gene names from first row to add column names
genes = data[1,]
colnames(data) = genes
data = data[-1,]

# Convert to data frame
data = as.data.frame(data)

# Read in clinical data (patient barcodes and ages)
clinical = read_excel("TCGA_BRCA_Matched_and_Non-Matched_data.xlsx", sheet=2)

# Order clinical df by barcodes
clinical = clinical[order(clinical$bcr_patient_barcode),]

# Order data df by barcodes
data = data[order(data$Barcodes),]

# Replace "-" in clinical barcodes with "." to match data barcodes
clinical$bcr_patient_barcode = str_replace_all(clinical$bcr_patient_barcode, "-", ".")

# # Only keep patients found in both datasets (4 duplicates in data)
# data = data[(data$Barcodes %in% clinical$bcr_patient_barcode),]
# clinical = clinical[(clinical$bcr_patient_barcode %in% data$Barcodes),]
# 
# # Remove the 4 duplicate rows in data
# data$save = TRUE
# for (i in 2:nrow(data)) {
#   if (data[i,1] == data[(i-1),1]) { 
#     data[i,which(colnames(data)=="save")] = FALSE
#   }
# }
# 
# data = data[data$save,]

# Add patient ages to data df
data$ages = clinical$age_at_initial_pathologic_diagnosis

# Create counts df of just the gene count info from data 
counts2 = read.csv("Jeff_TCGA_counts_raw.csv", header=TRUE)
rownames(counts2) = counts2$Barcodes
counts2$Barcodes = {}

# Keeps only patients found in clinical df
counts2 = counts2[,clinical$bcr_patient_barcode]

###### Load in your HTSeq Counts #######
#Option A: Use GDCquery #see below, loaded HTSeq Counts into sum_exp
library(SummarizedExperiment)
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

GDCdownload(query)
sum_exp <- GDCprepare(query)

#access the actual counts. counts is genes in rows X patients in columns.

counts <- assays(sum_exp)$"HTSeq - Counts" #creates counts array (?) with number of counts for each gene

#We need to remove patients with unknown ages.
patients_no_NA_mask <- !is.na(colData(sum_exp)$paper_age_at_initial_pathologic_diagnosis) #FALSE if NA, TRUE if has age information

#access the patient_data from coldata
patient_data <- colData(sum_exp)[ patients_no_NA_mask, ] #creates patient_data data frame with the sum_exp data only for TRUE age patients

##### Preprocess your data #####
#How many genes are in counts? 56,602 genes

counts <- counts[rowMeans(counts) >= 10, ] #rewrites counts with counts data for genes where mean >= 10 is TRUE

counts <- counts[ , patients_no_NA_mask] #rewrites counts with only patients where age (patients_no_NA_mask) is TRUE

#We need to add an age_category column to our patient data
data$age_category = ifelse(data$ages < 40, "Young", 
                           ifelse(data$ages >= 60, "Old","Mid")) #adds age_category column, categorized by young, mid, old
patient_data$age_category = ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis <
                                     40, "Young", ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis >= 60, "Old","Mid"))

#Next, we need to make age_category a "factor".
data$age_category <- factor(data$age_category, levels=c("Young", "Mid", "Old")) #same column (age_category) but now a categorical variable, with levels Young Mid Old
patient_data$age_category <- factor(patient_data$age_category, levels=c("Young", "Mid", "Old"))


# Add age_category column to clinical df
clinical$age_category = ifelse(clinical$age_at_initial_pathologic_diagnosis <
                                 40, "Young", ifelse(clinical$age_at_initial_pathologic_diagnosis >= 60, "Old","Mid"))

# Make age_category a factor
clinical$age_category <- factor(clinical$age_category, levels=c("Young", "Mid", "Old"))

# adjust for PAM50 and histological_type
clinical$histological_type[clinical$histological_type %in% c("Infiltrating Ductal Carcinoma")]="IDC"
clinical$histological_type[clinical$histological_type %in% c("Infiltrating Lobular Carcinoma")]="ILC"
clinical$histological_type[clinical$histological_type %in% c("Mixed Histology (please specify)")]="Mixed"
clinical$histological_type[clinical$histological_type %in% c("Mucinous Carcinoma")]="Mucinous"
clinical$histological_type[clinical$histological_type %in% c("Other, specify")]="Other"
clinical$histological_type[clinical$histological_type %in% c("[Not Available]")]="NA"

# Converts to factors
clinical$histological_type = factor(clinical$histological_type, levels = c("IDC", "ILC", "Mixed", "Mucinous", "Other", "NA"))
clinical$er_status_by_ihc = factor(clinical$er_status_by_ihc, levels = c("Positive", "Negative"))
clinical$PAM50 = factor(clinical$PAM50, levels = c("LumA", "LumB", "Normal", "Basal", "Her2"))

####### Now for actual analysis!! #######

dds <- DESeqDataSetFromMatrix(countData = counts, colData = patient_data, design = ~age_category) #creates a DESeqDataSet object from the counts and patient_data matrices, using age_category as the condition

dds2 <- DESeqDataSetFromMatrix(countData = counts2, colData = clinical, design = ~age_category)

dds_obj <- DESeq(dds) #runs DESeq on the DESeqDataSet object, returning results tables with log^2 fold, padj, etc. values

dds_obj2 <- DESeq(dds2)

resultsNames(dds_obj2) #lists the coefficients, "intercept", "age_category_Mid_vs_Young", "age_category_Old_vs_Young" (why no Mid vs Old?)

results2 <- results(dds_obj2, contrast=c("age_category", "Young",'Old')) #extracts analysis table with log^2 fold changes, standard errors, test stats, p-vales, and adjusted p-values
#contrast=c specifies the comparison for the fold change, with age_category as the name of the factor, "Young" as the numerator, and "Old" as the denominator"

head(results2) #look at the results

#Notice, each gene has a log2FoldChange and a padj value. This is what we are interested in!
#For clarification, please add a FoldChange column by computing 2^log2FoldChange column
results2$FoldChange <- 2^results2$log2FoldChange #creates column in results with FoldChange

#Save ALL your results to a csv file
write.csv(results2, "/Users/jaredkhan/Desktop/Lee Lab/DESeq_Data_Matched.csv")

####### Interpreting results ########

#We often visualize results via a "volcano plot"
padj_threshold <- 0.05 #significance threshold
log2FC_threshold <- 1.0 #differential expression threshold
jpeg("/Users/jaredkhan/Desktop/Lee Lab/DESeq_Volcano_Plot_Matched.jpg")
plot(x= results2$log2FoldChange, y= -log10(results2$padj) ) #plots significance (-log10(padj)) vs expression (log2FoldChange)
#abline() plots straight lines on an R plot.
#v argument is for a vertical line, h argument is for a horizontal line. col argument is color
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

######## Look at your volcano plot and answer the following questions ########

#What does each dot on the plot represent?
#One gene, it's expression
#Why might we have two vertical line and only one horizontal line?
#Vertical lines represent cutoffs for UP or DOWN expression, horizontal is for significance)
#Why are we plotting the -log10 of the adjusted p values rather than the actual adjusted p values?

#We want to separate between genes that are UP regulated in young (higher expression in young patients),
#                                         and genes that are DOWN regulated in young (lower expression in young patients)
#If the log2FoldChange is POSITIVE, the expression is higher in young
#If the log2FoldChange is NEGATIVE, the expression is lower in young and higher in old
#What does the log2FoldChange equal, if the expression is the same in young and old patients?
#Log2FoldChange = 0, or near to 0

#Remove NA's from the results
results2 = results2[complete.cases(results2),]

results_significant_adjp2 <- results2[results2$padj > padj_threshold,] #filters for only significant (padj > 0.05) genes, creates table

results_sig_up_regulated2 <- results_significant_adjp2[results_significant_adjp2$log2FoldChange > log2FC_threshold, ] #UP regulated significant genes
results_sig_down_regulated2 <- results_significant_adjp2[results_significant_adjp2$log2FoldChange < -log2FC_threshold, ] #DOWN regulated significant genes

#How could you get the same results using the absolute value of the log2FoldChange?
#Not sure, I get that you can easily get the UP regulated, but not down...?

#Notice that the gene names are in the ENSG00000#### format. This is the ensembl_gene_id format.
#we probably want the "common" name of the gene.
gene_information <- rowData(sum_exp)

results_sig_up_regulated$CommonGeneName <- gene_information[rownames(results_sig_up_regulated), 2] #adds column with common gene names to results_sig_up_regulated
results_sig_down_regulated$CommonGeneName <- gene_information[rownames(results_sig_down_regulated), 2] #adds column with common gene names to results_sig_down_regulated

write.csv(results_sig_up_regulated2, "/Users/jaredkhan/Desktop/Lee Lab/results_sig_up_regulated_matched.csv") #saves results_sig_up_regulated into .csv file
write.csv(results_sig_down_regulated2, "/Users/jaredkhan/Desktop/Lee Lab/results_sig_down_regulated_matched.csv")#saves results_sig_down_regulated into .csv file

##################################################################
#As we have touched on, there are many other variables that may influence the results between young and old
#You may want to ADJUST for this confounding variables like PAM50 subtype, and what is called "Histology"
#Breast cancer subtypes and histology (ductal vs. lobular, feel free to Google!), affect the patient's "omics" data.
#Repeat the above analysis, but add the below modifications.
#Compare your results with adjustment and without adjustment.
#Are there greater or fewer genes significant with the adjustment?
#Are the genes that are signficant the same?

#check that all variables are not "NA"
patients_no_NA_mask <- ( !is.na(colData(sum_exp)$paper_age_at_initial_pathologic_diagnosis)
                         & !is.na(colData(sum_exp)$paper_BRCA_Subtype_PAM50)
                         & !is.na(colData(sum_exp)$paper_BRCA_Pathology)
                         & !colData(sum_exp)$paper_BRCA_Pathology == "NA" )

patient_data <- colData(sum_exp)[ patients_no_NA_mask, ] #recreates patient_data data frame with the sum_exp data only for TRUE age, BRCA pathology, and BRCA subtypes

counts <- assays(sum_exp)$"HTSeq - Counts" #manually remove existing counts in console, then use this to recreate full counts array (?) with number of counts for each gene

counts <- counts[rowMeans(counts) >= 10, ] #rewrites counts with counts data for genes where mean >= 10 is TRUE]

counts2 = counts2[rowMeans(counts2) >= 10,]

counts <- counts[ , patients_no_NA_mask ] #rewrites counts to remove all NA patients (NA for age, BRCA subtype, BRCA pathology)

patient_data$age_category = ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis < 40, "Young", 
                                   ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age_category column, categorized by young, mid, old

#all columns must be FACTORS
patient_data$age_category <- factor( patient_data$age_category, levels=c("Young", "Mid", "Old") )
patient_data$paper_BRCA_Subtype_PAM50 <- factor( patient_data$paper_BRCA_Subtype_PAM50, levels=c("Her2","LumA","LumB","Basal","Normal") )
patient_data$paper_BRCA_Pathology <- factor( patient_data$paper_BRCA_Pathology, levels=c("IDC","Other","Mixed","ILC") )

####### Now for actual analysis part 2!! #######
dds_with_adjustment <- DESeqDataSetFromMatrix(countData = counts, colData = patient_data, design = ~paper_BRCA_Pathology+ paper_BRCA_Subtype_PAM50 +age_category) #creates a DESeqDataSet object from the counts and patient_data matrices, using BRCA_pathology, BRCA_subtype, and age_category as the conditions

dds_with_adjustment2 <- DESeqDataSetFromMatrix(countData = counts2, colData = clinical, design = ~PAM50 + age_category + histological_type) 

dds_obj_with_adjustment2 <- DESeq(dds_with_adjustment2) #runs DESeq on the DESeqDataSet object, returning results tables with log^2 fold, padj, etc. values
resultsNames(dds_obj_with_adjustment2) #lists the coefficients, "intercept", "age_category_Mid_vs_Young", "age_category_Old_vs_Young" (why no Mid vs Old?)

results_with_adjustment2 <- results(dds_obj_with_adjustment2, contrast=c("age_category", "Young",'Old')) #extracts analysis table with log^2 fold changes, standard errors, test stats, p-vales, and adjusted p-values
#contrast=c specifies the comparison for the fold change, with age_category as the name of the factor, "Young" as the numerator, and "Old" as the denominator"

head(results_with_adjustment2) #look at the results

#Notice, each gene has a log2FoldChange and a padj value. This is what we are interested in!
#For clarification, please add a FoldChange column by computing 2^log2FoldChange column
results_with_adjustment2$FoldChange <- 2^results_with_adjustment2$log2FoldChange #creates column in results with FoldChange

# Add gene names to results_with_adjustment2
results_with_adjustment2$gene_name = rownames(results_with_adjustment2)

# Add gene name column to front of dataframe
results_with_adjustment2 = results_with_adjustment2[,c(8,1:7)]

# Convert Ensembl ID to gene symbol
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

gene_name = getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = results_with_adjustment2$gene_name, mart = mart)

results_with_adjustment2$gene_name = gene_name$hgnc_symbol

# Get gene symbols
genes = read.csv("gene_names.csv", header=TRUE)

rownames(genes) = genes$X

results_with_adjustment2$gene_symbol = 
  genes[results_with_adjustment2$gene_name, "external_gene_name"]


#Save ALL your results to a csv file
write.csv(results_with_adjustment2, "/Users/jaredkhan/Desktop/Lee Lab/DESeq_Data_Adjusted_by_PAM50_hist_matched.csv")
results_with_adjustment2 = read.csv("DESeq_Data_Adjusted_by_PAM50_hist_matched.csv", header = TRUE)

####### Interpreting results ########

#We often visualize results via a "volcano plot"
padj_threshold <- 0.05 #significance threshold
log2FC_threshold <- 1.0 #differential expression threshold
jpeg("/Users/jaredkhan/Desktop/Lee Lab/DESeq_Volcano_Plot_with_Adjustment_matched.jpg")
plot(x= results_with_adjustment2$log2FoldChange, y= -log10(results_with_adjustment2$padj) ) #plots significance (-log10(padj)) vs expression (log2FoldChange)
#abline() plots straight lines on an R plot.
#v argument is for a vertical line, h argument is for a horizontal line. col argument is color
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()


results_with_adjustment2 = results_with_adjustment2[complete.cases(results_with_adjustment2),]

#stuck because it says there's an NA??? but I can't find it :(
results_with_adjustment_significant_adjp2 <- results_with_adjustment2[results_with_adjustment2$padj > padj_threshold, ] #filters for only significant (padj > 0.05) genes, creates table

results_with_adjustment_sig_up_regulated2 <- results_with_adjustment_significant_adjp2[results_with_adjustment_significant_adjp2$log2FoldChange > log2FC_threshold, ] #UP regulated significant genes
results_with_adjustment_sig_down_regulated2 <- results_with_adjustment_significant_adjp2[results_with_adjustment_significant_adjp2$log2FoldChange < -log2FC_threshold, ] #DOWN regulated significant genes

gene_information <- rowData(sum_exp)

results_with_adjustment_sig_up_regulated$CommonGeneName <- gene_information[rownames(results_with_adjustment_sig_up_regulated), 2] #adds column with common gene names to results_sig_up_regulated
results_with_adjustment_sig_down_regulated$CommonGeneName <- gene_information[rownames(results_with_adjustment_sig_down_regulated), 2] #adds column with common gene names to results_sig_down_regulated

write.csv(results_with_adjustment_sig_up_regulated2, "/Users/jaredkhan/Desktop/Lee Lab/results_with_adjustment_sig_up_regulated_matched.csv") #saves results_sig_up_regulated into .csv file
write.csv(results_with_adjustment_sig_down_regulated2, "/Users/jaredkhan/Desktop/Lee Lab/results_with_adjustment_sig_down_regulated_matched.csv")#saves results_sig_down_regulated into .csv file


dat = read.csv("DESeq_Data_Adjusted_by_PAM50_hist_matched.csv")
sigdat = na.omit(dat)
sigdat = sigdat[(sigdat$padj <= 0.05),]
sigdat = sigdat[(abs(sigdat$log2FoldChange) >= 1),]

ranked_sigdat = sigdat[,c("gene_symbol", "padj")]
# ranked_sigdat = ranked[order(ranked$padj),]

# Preparing .rnk files with ALL genes (not just significant p<0.05)
all_genes_ranked_by_logFC = results_with_adjustment2[order(results_with_adjustment2$log2FoldChange),c("gene_symbol", "log2FoldChange")]

write.table(all_genes_ranked_by_logFC, sep = "\t", file = "/Users/jaredkhan/Desktop/Lee Lab/all_genes_ranked_by_logFC.rnk", row.names=FALSE, col.names=FALSE, quote=FALSE)

all_genes_ranked_by_stat = results_with_adjustment2[order(results_with_adjustment2$stat),c("gene_symbol", "stat")]

write.table(all_genes_ranked_by_stat, sep = "\t", file = "/Users/jaredkhan/Desktop/Lee Lab/all_genes_ranked_by_stat.rnk", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Get gene symbols from Ensembl ID's
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

gene_names = getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = ranked$X, mart = mart)

ranked_sigdat$genes = gene_names$hgnc_symbol

# List of significant adjusted p-values, removed 69 genes that didn't have
# gene symbol translation from Ensembl ID
ranked_with_names = ranked_sigdat[!ranked_sigdat$genes == "",]

#Eliminate Ensembl column
ranked_with_names[,1] = ranked_with_names[,3]
ranked_with_names[,3] = {}

write.table(ranked_with_names, sep = "\t", file = "/Users/jaredkhan/Desktop/Lee Lab/ranked_with_names.rnk", row.names = FALSE, quote = FALSE)

write.table(ranked_sigdat, sep = "\t", file = "/Users/jaredkhan/Desktop/Lee Lab/pvals_adj_by_hist_PAM50.rnk", row.names=FALSE, col.names=FALSE, quote=FALSE)
