########################################################################
########################################################################
# Prepare table 1 summary for each dataset in the analysis.
# Requires metagxbreast package to be in the present directory;
# selects studies using the GBJ results in reference files folder (also in pwd).
########################################################################
########################################################################

library(xtable)

GBJ_merged_results <- read.delim("Reference_files/GBJ_merged_results.txt", stringsAsFactors=FALSE)
studies <- names(GBJ_merged_results)[3:23]

outMatrix <- matrix(NA, nrow=21, ncol=6)
j <-0
for (study in studies){
j <- j+1
  if (study == "GSE25066"){
  	#Assume metagxbreast package is in same directory as reference files folder
    load(paste0("MetaGxBreast/data/",study,".RData"))
  } else {
    load(paste0("MetaGxBreast/data/",study,".rda"))
  }
 phenData <- (get(study)@phenoData@data)
 phenotypic <- phenData[phenData$grade %in% c(1,3),]
 phenotypic$gradeInd <- ifelse(phenotypic$grade==3, 1, 0)
 outMatrix[j, 1] <- study
 outMatrix[j, 2] <- dim(phenotypic)[1]
 outMatrix[j, 3] <- paste(sum(phenotypic$gradeInd), ", (", 100*round(sum(phenotypic$gradeInd)/length(phenotypic$gradeInd),3), "%)", sep="")
 outMatrix[j, 4] <- paste(sum(phenotypic$er=="positive", na.rm = T), ", (", 100*round(sum(phenotypic$er=="positive", na.rm = T)/sum(is.na(phenotypic$er)==FALSE),3), "%)", sep="")
 outMatrix[j, 5] <- paste(sum(phenotypic$her2=="positive", na.rm = T), ", (", 100*round(sum(phenotypic$her2=="positive", na.rm = T)/sum(is.na(phenotypic$her2)==FALSE),3), "%)", sep="")
 outMatrix[j, 6] <- round(mean(phenotypic$age_at_initial_pathologic_diagnosis, na.rm=TRUE),1)
 }
xtable(outMatrix)
