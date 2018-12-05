########################################################################
########################################################################
# Run GSA for each dataset in the compendium
# Libraries required
library(GSA);library(stringr)
########################################################################
########################################################################

#Read in the gmt file of gene sets (V5.2 TF set from msigdb)
geneSetFull <- read.csv('GSEA_data/tfSets.csv')
geneset.names <- as.character(geneSetFull[,2])
genesetsBuild <- geneSetFull[, 3:2942]
#Build list of gene sets 
genesets <- vector("list",nrow(genesetsBuild))
for(i in 1:nrow(genesetsBuild)){
  genesets[[i]]=paste("geneid.",as.character(genesetsBuild[i,])[as.character(genesetsBuild[i,])!="NA"],sep="")
}

####################################
#### Run for each study
####################################
studyFiles <- list.files("GSEA_data/")
studies <- gsub( "_.*$", "", studyFiles[which(str_sub(studyFiles, start= -3)=="cls")] )
for(aID in studies){
  #Read in the grade information for outcome y
  pheno <- read.delim(paste0("GSEA_data/",aID,"_pheno.cls"),skip=1,sep=" ")
  y <- as.numeric(unname(pheno[1,]))+1
  
  #Read in the observations
  x <- as.matrix(read.delim(paste0("GSEA_data/",aID,"_exp.txt"), row.names=1, stringsAsFactors=FALSE))
  genenames <- row.names(x)
  
  #Run GSA for all of the sets
  GSA.obj<-GSA(x,y, genenames=genenames, genesets=genesets, resp.type="Two class unpaired", nperms=5000)
  GSA.out <- rbind(GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=1)$negative, 
                   GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=1)$positive)
  write.table(GSA.out, paste0('GSA_out/GSA_results_',aID,'.txt'), append=F, quote=F, row.names=F, col.names=T, sep='\t') 
}



#############################################################################
#############################################################################
#Organize GSA data for output

# This directory holds all GSA results (and nothing else)
GSA_files_dir <- 'RevisionCode/GSA'
ref_dir <- 'Reference_files/'


setwd(GSA_files_dir)
fname_vec <- list.files()
period_pos <- regexpr('[.]', fname_vec)
fname_roots <- substr(fname_vec, 13, period_pos-1)

# Loop through files and record test stat/repeat for p-value
results <- NULL
case_control_tab <- data.frame(Study=fname_roots, Cases=NA, Controls=NA)
for (i in 1:length(fname_vec)) {
  
  # Read the file
  temp_tab <- tryCatch(read.table(fname_vec[i], sep='\t', header=F), warning=function(w) w, error=function(e) e)
  if (length(class(temp_tab)) > 1) {next}
  temp_tab <- temp_tab[2:dim(temp_tab)[1],]
  
  # Parse file name, record name
  underscore_pos <- gregexpr('_', fname_vec[i])[[1]]
  num_cases <- as.numeric(substr(fname_vec[i], underscore_pos[3]+1, underscore_pos[4]-1))
  num_controls <- as.numeric(substr(fname_vec[i], underscore_pos[5]+1, underscore_pos[6]-1))
  case_control_tab$Cases[i] <- num_cases
  case_control_tab$Controls[i] <- num_controls
  
  # First one
  if (is.null(results)) {
    results <- temp_tab
    results <- results[,c(2,3)]
    colnames(results) <- c('V2', fname_roots[i])
  } else {
    results <- merge(results, temp_tab[,2:3], by="V2", all=T)
    colnames(results)[ncol(results)] <- fname_roots[i]
  }
}
colnames(results)[1] <- 'pathway'


##################################################################################
# Add pathway size to results
setwd(ref_dir)
pathways_tested <- read.table('Breast_tf_pathways_120116.txt')
pathways_with_size <- read.table('Breast_tf_pathways_120116_size.txt', header=T)
merged_results <- merge(results, pathways_with_size, by='pathway')
merged_results <- merged_results[, c(1,ncol(merged_results), 2:(ncol(merged_results)-1))]
nrow(merged_results)

# Write it
write.table(merged_results, 'GSA_merged_results_stat.txt', append=F, quote=F, row.names=F, col.names=T, sep='\t')


# This directory holds all GSA results (and nothing else)
GSA_files_dir <- 'RevisionCode/GSA'
ref_dir <- 'Reference_files/'


setwd(GSA_files_dir)
fname_vec <- list.files()
period_pos <- regexpr('[.]', fname_vec)
fname_roots <- substr(fname_vec, 13, period_pos-1)

# Loop through files and record test stat/repeat for p-value
results <- NULL
case_control_tab <- data.frame(Study=fname_roots, Cases=NA, Controls=NA)
for (i in 1:length(fname_vec)) {
  
  # Read the file
  temp_tab <- tryCatch(read.table(fname_vec[i], sep='\t', header=F), warning=function(w) w, error=function(e) e)
  if (length(class(temp_tab)) > 1) {next}
  temp_tab <- temp_tab[2:dim(temp_tab)[1],]
  
  # Parse file name, record name
  underscore_pos <- gregexpr('_', fname_vec[i])[[1]]
  num_cases <- as.numeric(substr(fname_vec[i], underscore_pos[3]+1, underscore_pos[4]-1))
  num_controls <- as.numeric(substr(fname_vec[i], underscore_pos[5]+1, underscore_pos[6]-1))
  case_control_tab$Cases[i] <- num_cases
  case_control_tab$Controls[i] <- num_controls
  
  # First one
  if (is.null(results)) {
    results <- temp_tab
    results <- results[,c(2,4)]
    colnames(results) <- c('V2', fname_roots[i])
  } else {
    results <- merge(results, temp_tab[,2:4], by="V2", all=T)
    colnames(results)[ncol(results)] <- fname_roots[i]
  }
}
colnames(results)[1] <- 'pathway'


##################################################################################
# Add pathway size to results
setwd(ref_dir)
pathways_tested <- read.table('Breast_tf_pathways_120116.txt')
pathways_with_size <- read.table('Breast_tf_pathways_120116_size.txt', header=T)
merged_results <- merge(results, pathways_with_size, by='pathway')
merged_results <- merged_results[, c(1,ncol(merged_results), 2:(ncol(merged_results)-1))]
nrow(merged_results)

# Write it
write.table(merged_results, 'GSA_merged_results_pval.txt', append=F, quote=F, row.names=F, col.names=T, sep='\t')



