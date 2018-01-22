########################################################################
########################################################################
# This code prepares the MetaGxBreast data so that it is in a format
# that can be analyzed by GSEA

# ER status positive or negative
desired_er <- 'positive'
# Where do you want to save the GSEA-formatted files?
GSEA_data_dir <- '/users/ryansun/desktop/GSEA_data'

########################################################################
########################################################################


########################################################################
# Data resides in this package
library(MetaGxBreast)

# Load in the cleaned expression data
source(system.file("extdata", "patientselection.config", package="MetaGxBreast"))
source(system.file("extdata", "createEsetList.R", package="MetaGxBreast"))

# Write the data
for (dataset_it in 1:length(names(esets))) {
	temp_eset_name <- names(esets)[dataset_it]
	expData <- as.data.frame(esets[[dataset_it]])
	
	# These have lots of NA
	if (temp_eset_name %in% c('NCI', 'NKI', 'UCSF')) {next}

	# Keep only the correct ER status	
	expData <- expData[which(expData$er == desired_er), ]
	
	# Only pick out grades 1 and 3
	expData <- expData[which(expData$grade %in% c(1,3)),]
	expData$gradeInd <- ifelse(expData$grade==3, 1, 0)
	
	# Too few subjects qualifying? Next
	if (nrow(expData) < 10) {next}
	num_cases <- length(which(expData$gradeInd == 1))
	num_controls <- length(which(expData$gradeInd == 0))
	if (num_cases < 5 | num_controls < 5) {next}
	
	# Figure out which columns in expData correspond to genes
	temp_colnames <- substr(names(expData), 1, 7)
	gene_cols <- which(temp_colnames == 'geneid.')
	
	# This is the data we will save for GSEA
	data_to_write <- t(expData[, gene_cols])
	
	# Have to make the gene names a column, so we can give them a column heading
	data_to_write <- data.frame(NAME=rownames(data_to_write), data_to_write)
	
	# If there are NAs, will have problems
	any_NA <- which(is.na(data_to_write)) 
	bad_cols <- rep(NA, ncol(data_to_write))
	if (length(any_NA) > 0) {
		for (col_it in 2:ncol(data_to_write)) {
			temp_col <- data_to_write[, col_it]
			if (length(which(is.na(temp_col))) > 0) {
				bad_cols[col_it] <- 1
			}
		}
		data_to_write <- data_to_write[, -which(bad_cols == 1)]
		
		# Tell us how many we had to remove
		cat(temp_eset_name, "removed ", length(bad_cols), " subjects for NA and ", ncol(data_to_write), " left \n")
	}

	# Save it
	setwd(GSEA_data_dir)
	temp_fname <- paste(temp_eset_name, '_exp.txt', sep='')
	write.table(data_to_write, file=temp_fname, sep='\t', quote=F, append=F, row.names=F, col.names=T)	
	
	# Make phenotype file
	pheno_row3 <- expData$gradeInd
	if (length(any_NA) > 0) {
		pheno_row3 <- pheno_row3[-which(bad_cols==1)]
	}
	pheno_row1 <- c(length(pheno_row3), 2, 1, rep(NA, length(pheno_row3)-3))
	if (pheno_row3[1] == 1) {
		pheno_row2 <- c('#', 'Case', 'Control', rep(NA, length(pheno_row3)-3))
	} else {
		pheno_row2 <- c('#', 'Control', 'Case', rep(NA, length(pheno_row3)-3))
	}
	
	# Write it
	pheno_file <- rbind(pheno_row1, pheno_row2, pheno_row3)
	temp_fname <- paste(temp_eset_name, '_pheno.cls', sep='')
	write.table(pheno_file, file=temp_fname, sep=' ', append=F, quote=F, row.names=F, col.names=F, na='')
	
}


############################################################################
############################################################################
# At this point, you need to use one of the GSEA utilities to run the GSEA
# analysis on all of the .txt and .cls files.  The .gmt file is given
# by msigdb_tf.gmt
############################################################################
############################################################################


output_dir <- '/users/ryansun/desktop/plot_data'

########################################################################
########################################################################
# GSEA will generally output all your results into a folder
# that is labeled with the day's date, e.g. /sep02/.
# Give that filepath here (with the backslash at the end).
GSEA_root <- '/users/ryansun/gsea_home/output/sep02/'

########################################################################
########################################################################


########################################################################
# We need to get the pertinent information from the huge output of GSEA
# and save just what we need.
library(XML)
setwd(GSEA_root)

# Get all results directory names
all_folders <- list.files()

# Loop through folders, should have one study results in each folder
for (folder_it in 1:length(all_folders)) {
	
	cat("Starting ", folder_it, '\n')
	
	temp_append <- all_folders[folder_it]
	temp_edb <- paste(GSEA_root, temp_append, '/edb', sep='')
	temp_dir <- paste(GSEA_root, temp_append, sep='')
	setwd(temp_edb)
	
	# Figure out which study
	edb_names <- list.files()
	edb_ext <- substr(edb_names, nchar(edb_names)-2, nchar(edb_names))
	cls_file <- which(edb_ext == 'cls')
	cls_name <- edb_names[cls_file]
	# Figure out where the underscore is
	underscore_pos <- regexpr('_', cls_name)[1]
	study_name <- substr(cls_name, 1, underscore_pos-1)
	
	# Get the names of the reports
	period_pos <- gregexpr('\\.', temp_append)
	study_num <- substr(temp_append, period_pos[[1]][2]+1, nchar(temp_append))
	
	# Name the two study results tables
	tab0_name <- paste('gsea_report_for_0_', study_num, '.html', sep='')
	tab1_name <- paste('gsea_report_for_1_', study_num, '.html', sep='')
	
	# Open the tables
	setwd(temp_dir)
	tab0 <- readHTMLTable(tab0_name, header=F)[[1]]
	tab1 <- readHTMLTable(tab1_name, header=F)[[1]]
	
	# Column 6 is NES and 9 is FWER p-val
	# Put them together, sort by FWER p-val
	tot_results <- rbind(tab0, tab1)
	tot_results <- tot_results[order(tot_results$V9, decreasing=FALSE), ]
	
	# Write to GSEA data folder
	tot_results <- tot_results[, c(2,4:9)]
	colnames(tot_results) <- c('Pathway', 'Num_genes', 'ES', 'NES', 'Nom_p', 'FDR', 'FWER')
	setwd(GSEA_data_dir)
	fname_out <- paste(study_name, '_GSEA_results.txt', sep='')
	write.table(tot_results, fname_out, append=F, quote=F, row.names=F, col.names=T, sep='\t')
}
########################################################################
########################################################################



########################################################################
# Merge all the GSEA results into one file.
setwd(GSEA_data_dir)

all_GSEA_studies <- c('CAL', 'DFHCC', 'DFHCC2', 'EXPO', 'GSE25066', 'GSE32646',
					'GSE58644', 'IRB', 'MAINZ', 'MAQC2', 'METABRIC', 'PNC', 'STK',
					'STNO2', 'TRANSBIG', 'UNC4', 'UNT', 'UPP')
all_GSEA_fnames <- paste(all_GSEA_studies, '_GSEA_results.txt', sep='')

GSEA_results <- NULL
for (i in 1:length(all_GSEA_fnames)) {
	
	temp_file <- read.table(all_GSEA_fnames[i], header=T)
	
	if (is.null(GSEA_results)) {
		GSEA_results <- temp_file[, c(1,2,7)]
		colnames(GSEA_results)[ncol(GSEA_results)] <- all_GSEA_studies[i]
	} else {
		temp_file <- temp_file[, c(1, ncol(temp_file))]
		GSEA_results <- merge(GSEA_results, temp_file, by='Pathway', all.x=T)
		colnames(GSEA_results)[ncol(GSEA_results)] <- all_GSEA_studies[i]
	}
}

# Record merged results
setwd(GSEA_data_dir)
write.table(GSEA_results, 'GSEA_merged_results.txt', append=F, quote=F, row.names=F, col.names=F, sep='\t')
