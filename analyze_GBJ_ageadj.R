########################################################################
########################################################################
# First prepare the data for analysis
# Seeking to run GBJ with age as covariate

########################################################################
########################################################################
studyFiles <- list.files("~/Desktop/GSEA_data/")
studies <- gsub( "_.*$", "", studyFiles[which(str_sub(studyFiles, start= -3)=="cls")] )

for (study in studies){
  if (study == "GSE25066"){
  	#Location of metagxbreast package
    load(paste0("~/Desktop/MetaGxBreast/data/",study,".RData"))
  } else {
    load(paste0("~/Desktop/MetaGxBreast/data/",study,".rda"))
  }
  #Read in relevant exp/phen
  expData <- as.data.frame(t(read.delim(paste0("~/Desktop/GSEA_data/",study,"_exp.txt"), row.names=1, stringsAsFactors=FALSE)))
  phenData <- as.matrix(get(study)@phenoData@data)
  #Get the observations who have expression we use
  phenData_Match <- data.frame(phenData[which(rownames(phenData) %in% rownames(expData)),])
  #double check the names line up 
  table(rownames(phenData_Match) == rownames(expData))
  #write out age file
  write.table(phenData_Match$age_at_initial_pathologic_diagnosis, paste0('~/Desktop/SensitivityData/',study,'_age.txt'), append=F, quote=F, row.names=T, col.names=F, sep='\t') 
}


#Get additional studies age info
studies <- c("NCI", "NKI", "UCSF")
for (study in studies){
  load(paste0("~/Desktop/MetaGxBreast/data/",study,".rda"))
  table(is.na(phenoData(get(study))$age[ (phenoData(get(study))$er) == 'positive' & (phenoData(get(study))$grade %in% c(1,3)) ]))
}


########################################################################
########################################################################
# Run GBJ for each dataset in the compendium (run on cluster)
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
ageData <- read.delim(paste0(currentData,"_age.txt"),header=F)
table(row.names(expData) == ageData$V1) #check matching obs
expData$age <- ageData$V2
expData <- expData[is.na(expData$age)==FALSE,] #remove those missing


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
null_mod <- glm(expData$gradeInd~expData$age, family=binomial(link = "logit"))

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

      #Run GBJ to get results
      null_mod <- glm(expData$gradeInd~expData$age, family=binomial(link = "logit"))
      reg_stats <- calc_score_stats(null_model=null_mod, factor_matrix=genotype_data, link_function="logit")
      old_test_stats <- reg_stats$test_stats
      correct_pvalues <- pnorm(-abs(old_test_stats), mean=emp_mu, sd=emp_sd_ns) +
                          (1-pnorm(abs(old_test_stats), mean=emp_mu, sd=emp_sd_ns))
      new_test_stats <- qnorm(correct_pvalues/2)
      GBJOut <- GBJ(test_stats=new_test_stats, cor_mat=reg_stats$cor_mat)
      return( c(as.character(geneSets[k,1]), GBJOut$GBJ, GBJOut$GBJ_pvalue, GBJOut$err_code) )
    }
  }

  return( c(as.character(geneSets[k,1]), NA, NA, NA) )
}

fname_out <- paste(currentData, 'cases', num_cases,
					'controls', num_controls, 'GBJ_age.txt', sep='_')
write.table(results, fname_out, append=F, quote=F, row.names=F, col.names=F, sep='\t')
