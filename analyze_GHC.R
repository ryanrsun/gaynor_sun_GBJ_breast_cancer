########################################################################
########################################################################
# Run GHC for each dataset in the compendium
# Takes one trailing argument (aID), which controls the dataset to analyze,
# you can also manually set this value.

# Libraries location
libs <- '/n/home13/rsun/Rlibrary/3.0.1'


# Read in aID
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
########################################################################
########################################################################

# Load necessary libraries
library(splines)
library(foreach)
library(doMC, lib= libs)
library(logging, lib= libs)
library(GBJ, lib= libs)

# Register for parallel
registerDoMC(4)

# Read the gene sets, get rid of first column (just numbering)
geneSets <- read.csv('tfSets.csv')
geneSets <- geneSets[, 2:2942]

# Load in the cleaned expression data
library(stringr)
studyFiles <- list.files()
studies <- gsub( "_.*$", "", studyFiles[which(str_sub(studyFiles, start= -3)=="cls")] )

# Pick out an expression data set by aID
currentData <- studies[aID]
print(aID)
expData <- as.data.frame(t(read.delim(paste0(currentData,"_exp.txt"), row.names=1, stringsAsFactors=FALSE)))
pheno <- read.delim(paste0(currentData,"_pheno.cls"),skip=1,sep=" ")
expData$gradeInd <- as.numeric(unname(pheno[1,]))

# Too few subjects qualifying? Next
if (nrow(expData) < 10) {stop()}
num_cases <- length(which(expData$gradeInd == 1))
num_controls <- length(which(expData$gradeInd == 0))
if (num_cases < 5 | num_controls < 5) {stop()}

# Before running test, get the empirical null distribution
# Need all the individual gene test statistics
all_genes <- unlist(geneSets[, 2:ncol(geneSets)])
all_genes <- unique(all_genes)
all_genes <- paste("geneid.", all_genes, sep="")
genotype_data <- expData[,which( names(expData) %in% all_genes )]

# Will also need this null model for later
null_mod <- glm(expData$gradeInd~1, family=binomial(link = "logit"))

# Just the score statistics for now to get empirical null
test_stats <- score_stats_only(null_model=null_mod, factor_matrix=genotype_data, link_function="logit")

# Remove any NAs
NA_ind <- which(is.na(test_stats))
if (length(NA_ind) > 0) {
  test_stats <- test_stats[-NA_ind]
}

# Find the endpoints
lower_pt <- round(min(test_stats), digits=1) - 0.1
upper_pt <- round(max(test_stats), digits=1) + 0.1

# Define the histogram bins, get the counts
hist_breaks <- seq(from=lower_pt, to=upper_pt, by=0.1)
counts <- hist(test_stats, breaks=hist_breaks)$counts
x_vec <- hist(test_stats, breaks=hist_breaks)$mids

# Smoothing matrix
M_mat <- matrix(data=NA, nrow=length(counts), ncol=length(counts))
lambda <- 1
for (i in 1:nrow(M_mat)) {
  temp_row <- dnorm( (rep(x_vec[i], length(x_vec)) - x_vec) / lambda ) / lambda
  c_k <- 1 / sum(temp_row)
  M_mat[i, ] <- temp_row * c_k
}

# Fit the model
carrier_vec <- M_mat %*% counts
ns_mod <- glm(counts ~ ns(x_vec, df=7),
              family=poisson(link='log'), offset=log(carrier_vec))

# Find mean
max_index <- which(ns_mod$fitted.values == max(ns_mod$fitted.values))
max_x <- x_vec[max_index]
emp_mu <- max_x

# Get the data for sd model
sd_ind <- which(x_vec >= max_x-1.51 & x_vec <= max_x+1.51)
sd_dat <- cbind(x_vec[sd_ind],
                ns_mod$fitted.values[sd_ind])
new_x <- sd_dat[,1]
new_xsquared <- sd_dat[,1]^2

# Natural spline 7 df
sd_mod_ns <- lm(log(sd_dat[,2] / (length(test_stats)*0.1)) ~ new_x + new_xsquared)
emp_sd_ns <- 1 / sqrt(-2*summary(sd_mod_ns)$coefficients[3,1])

# Now we can start running the actual analysis
results <- foreach(k=1:nrow(geneSets), .combine=rbind) %dopar% {
  genes <- geneSets[k,]
  genes <- genes[!is.na(genes)]
  genes <- paste("geneid.", genes[2:length(genes)], sep="")
  if (length(which( names(expData) %in% genes ))>1) {
    genotype_data <- expData[,which( names(expData) %in% genes )]

    if (dim(genotype_data)[2] > 1){

      # Check for NA
      column_sums <- colSums(genotype_data)
      if (length(which(is.na(column_sums))) > 0) {
        genotype_data <- genotype_data[, -which(is.na(column_sums))]
		
		# If only one column, ncol will be null
		if (is.null(ncol(genotype_data))) {
		  return( c(as.character(geneSets[k,1]), NA, NA, NA) )
		}
        if (ncol(genotype_data) < 2) {
          return( c(as.character(geneSets[k,1]), NA, NA, NA) )
        }
      }

      #Run GHC to get results
      null_mod <- glm(expData$gradeInd~1, family=binomial(link = "logit"))
      reg_stats <- calc_score_stats(null_model=null_mod, factor_matrix=genotype_data, link_function="logit")
      old_test_stats <- reg_stats$test_stats
      correct_pvalues <- pnorm(-abs(old_test_stats), mean=emp_mu, sd=emp_sd_ns) +
                          (1-pnorm(abs(old_test_stats), mean=emp_mu, sd=emp_sd_ns))
      new_test_stats <- qnorm(correct_pvalues/2)
      GHCOut <- GHC(test_stats=new_test_stats, cor_mat=reg_stats$cor_mat)
      return( c(as.character(geneSets[k,1]), GHCOut$GHC, GHCOut$GHC_pvalue, GHCOut$err_code) )
    }
  }

  return( c(as.character(geneSets[k,1]), NA, NA, NA) )
}

fname_out <- paste(currentData, 'cases', num_cases,
					'controls', num_controls, 'GHC.txt', sep='_')
write.table(results, fname_out, append=F, quote=F, row.names=F, col.names=F, sep='\t')








#############################################################################
#############################################################################
#Organize GHC data for output

# This directory holds all GHC results (and nothing else)
GHC_files_dir <- 'RevisionCode/GHC'
ref_dir <- 'Reference_files/'


setwd(GHC_files_dir)
fname_vec <- list.files()
period_pos <- regexpr('_', fname_vec)
fname_roots <- substr(fname_vec, 1, period_pos-1)

# Loop through files and record p-values
results <- NULL
case_control_tab <- data.frame(Study=fname_roots, Cases=NA, Controls=NA)
for (i in 1:length(fname_vec)) {
  
  # Read the file
  temp_tab <- tryCatch(read.table(fname_vec[i], sep='\t', header=F), warning=function(w) w, error=function(e) e)
  if (length(class(temp_tab)) > 1) {next}
  
  # Fix the <1*10^(-12)
  err_rows <- which(temp_tab[,4] != 0 & temp_tab[,2] > 25)
  if (length(err_rows) > 0) {
    cat(length(err_rows), 'very significant pathways in study', i, 'truncated to p=10^-12', '\n')
    #temp_tab[err_rows, 3] <- 1*10^(-12)
  }
  
  # Parse file name, record name
  underscore_pos <- gregexpr('_', fname_vec[i])[[1]]
  num_cases <- as.numeric(substr(fname_vec[i], underscore_pos[3]+1, underscore_pos[4]-1))
  num_controls <- as.numeric(substr(fname_vec[i], underscore_pos[5]+1, underscore_pos[6]-1))
  case_control_tab$Cases[i] <- num_cases
  case_control_tab$Controls[i] <- num_controls
  
  # First one
  if (is.null(results)) {
    results <- temp_tab
    results <- results[,c(1,3)]
    colnames(results) <- c('pathway', fname_roots[i])
  } else {
    results <- cbind(results, temp_tab[,3])
    colnames(results)[ncol(results)] <- fname_roots[i]
  }
}


##################################################################################
# Add pathway size to results
setwd(ref_dir)
pathways_tested <- read.table('Breast_tf_pathways_120116.txt')
pathways_with_size <- read.table('Breast_tf_pathways_120116_size.txt', header=T)
merged_results <- merge(results, pathways_with_size, by='pathway')
merged_results <- merged_results[, c(1,ncol(merged_results), 2:(ncol(merged_results)-1))]
nrow(merged_results)

# Write it
write.table(merged_results, 'GHC_merged_results.txt', append=F, quote=F, row.names=F, col.names=T, sep='\t')