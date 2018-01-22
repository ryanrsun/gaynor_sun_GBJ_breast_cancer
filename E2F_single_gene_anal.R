################################################################################
################################################################################
# Get the single-gene statistics for E2F transcription factors
# Produces E2F_results.txt

# Change output directory
output_dir <- "/Users/ryansun/Dropbox/Research/Paper5/Rep22-Dec6/Final_code"

################################################################################
################################################################################

# Load libraries
library(MetaGxBreast)
library(GBJ)
# Load in the cleaned expression data
source(system.file("extdata", "patientselection.config", package="MetaGxBreast"))
source(system.file("extdata", "createEsetList.R", package="MetaGxBreast"))

# The E2F family
E2F_genes <- c(1870, 1871, 1874:1876, 144455, 79733)
E2F_genes <- paste('geneid.', E2F_genes, sep='')
results <- matrix(data=NA, nrow=length(E2F_genes)*length(esets), ncol=4)
results <- data.frame(results)
colnames(results) <- c('Study', 'Gene', 'Pvalue', 'Percentile')

# Loop through datasets
record_row <- 1
for (dataset_it in 1:length(names(esets))) {
	cat("Starting ", dataset_it, '\n')
	dataset_name <- names(esets)[dataset_it]
	
	temp_eset_name <- names(esets)[dataset_it]
	expData <- as.data.frame(esets[[dataset_it]])
	
	# These have lots of NA
	# if (temp_eset_name %in% c('NCI', 'NKI', 'UCSF')) {next}

	# Only keep ER positive/negative
	desired_er <- 'positive'
	expData <- expData[which(expData$er == desired_er), ]
	
	# Only pick out grades 1 and 3
	expData <- expData[which(expData$grade %in% c(1,3)),]
	expData$gradeInd <- ifelse(expData$grade==3, 1, 0)
	
	# Too few subjects qualifying? Next
	if (nrow(expData) < 10) {next}
	num_cases <- length(which(expData$gradeInd == 1))
	num_controls <- length(which(expData$gradeInd == 0))
	if (num_cases < 5 | num_controls < 5) {next}
	
	# Pick out the E2F genes only
	E2F_cols <- which(names(expData) %in% E2F_genes)
	if (length(E2F_cols) == 0) {next}
	E2F_data <- as.matrix(expData[, E2F_cols])
	if (ncol(E2F_data) == 1) {
		temp_E2F_names <- names(expData[E2F_cols])
	} else {
		temp_E2F_names <- colnames(E2F_data)
	}
	
	
	# Outcome
	outcome <- expData$gradeInd
	
	# Remove NAs
	temp_rowsums <- rowSums(E2F_data)
	bad_sub <- which(is.na(temp_rowsums))
	if (length(bad_sub) > 0) {
		E2F_data <- E2F_data[-bad_sub, ]
		outcome <- outcome[-bad_sub]
	}
	
	# Calculate score statistics
	null_mod <- glm(outcome~1, family=binomial(link = "logit"))
	E2F_stats <- score_stats_only(null_model=null_mod, factor_matrix=E2F_data, link_function="logit")
	
	# Figure out which columns in expData correspond to genes
	expData_names <- names(expData)
	expData_names <- substr(expData_names, 1, 7)
	gene_columns <- which(expData_names == 'geneid.')
	all_gene_data <- expData[, gene_columns]
	
	# Get all the score statistics
	# New null mod, different way of removing NAs
	null_mod <- glm(expData$gradeInd~1, family=binomial(link = "logit"))
	all_stats <- score_stats_only(null_model=null_mod, factor_matrix=all_gene_data, link_function="logit")
	
	# Remove any NAs
	NA_ind <- which(is.na(all_stats))
	if (length(NA_ind) > 0) {
	  all_stats <- all_stats[-NA_ind]
	}
	all_pvalues <- 1-pchisq(all_stats^2, df=1)
	
	# Record
	for (j in 1:length(temp_E2F_names)) {
		results$Gene[record_row] <- temp_E2F_names[j]
		results$Pvalue[record_row] <- 1-pchisq(E2F_stats[j]^2, df=1)
		results$Percentile[record_row] <- length(which(all_pvalues < results$Pvalue[record_row])) / length(all_pvalues)
		results$Study[record_row] <- dataset_name
		
		record_row <- record_row + 1
	}
}	

# Write
results <- results[1:(record_row-1),]
setwd(output_dir)
write.table(results, 'E2F_results.txt', append=F, quote=F, row.names=F, col.names=T, sep='\t')
