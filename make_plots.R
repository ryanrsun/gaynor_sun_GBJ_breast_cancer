#############################################################################
#############################################################################
# Will output the heatmap, files necessary for EnrichmentMap network,
# and Table 2 for breast cancer gene sets paper.

# This directory holds all 22 single-study GBJ results (and nothing else)
GBJ_files_dir <- '/users/ryansun/desktop/gsea_data/er_positive/'

# This directory holds the reference files:
# Breast_tf_pathways_120116.txt
# Breast_tf_pathways_120116_size.txt
# case_control_counts.txt
# GSEA_merged_results.txt
ref_dir <- '/users/ryansun/desktop/gsea_data/ref_dir/'

# This is the output directory
output_dir <- '/users/ryansun/desktop/gsea_data/ref_dir/'

#############################################################################
#############################################################################



# There should only be the 22 results files in this directory. 
setwd(GBJ_files_dir)
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
setwd(output_dir)
write.table(merged_results, 'GBJ_merged_results.txt', append=F, quote=F, row.names=F, col.names=T, sep='\t')
##################################################################################



##################################################################################
# Perform a meta-analysis with Fisher's Method
fisher_meta <- data.frame(pathway=merged_results$pathway, num_genes=merged_results$num_genes, Fisher_stat=NA, Fisher_p=NA)
for (i in 1:nrow(merged_results)) {
	# Checkpointing
	if (i%%100 == 0) {cat(i)}

	# Remove studies not contributing
	temp_row <- merged_results[i, 3:ncol(merged_results)]
	bad_ind <- which(is.na(temp_row))
	if (length(bad_ind) >= (length(fname_vec) - 1)) {next}
	if (length(bad_ind) > 0) {temp_row <- temp_row[-bad_ind]}
	
	# Fisher statistic
	fisher_stat <- -2*sum(log(temp_row))
	fisher_p <- 1-pchisq(q=fisher_stat, df=2*length(temp_row))
		
	# Record
	fisher_meta$Fisher_stat[i] <- fisher_stat
	fisher_meta$Fisher_p[i] <- fisher_p
}
fisher_meta <- fisher_meta[order(fisher_meta$Fisher_p), ]
if (length(which(is.na(fisher_meta$Fisher_p))) > 0) {
	fisher_meta <- fisher_meta[-which(is.na(fisher_meta$Fisher_p)), ]
}
##################################################################################


# Write
setwd(output_dir)
write.table(fisher_meta, file='ERpositive_meta_fisher.txt', append=F, quote=F, row.names=F, col.names=T)
###########################################################


###########################################################			
# Parse meta analysis results for EnrichmentMap - Using Fisher results!
# The GMT file should have pathway_name <TAB> description <TAB> gene1 <TAB> gene2 <TAB> ...
meta_top_pathways <- fisher_meta[which(fisher_meta$Fisher_p < (0.05 / nrow(fisher_meta))), ]
EM_rows <- which(as.character(pathways_tested[,1]) %in% as.character(meta_top_pathways$pathway))
EM_pathways <- pathways_tested[EM_rows, ]

# Get the correct ordering of pathways
order_EM_pathways <- match(meta_top_pathways$pathway, EM_pathways[,1])
EM_pathways <- EM_pathways[order_EM_pathways, ]

# The results file should have pathway_name <TAB> description <TAB> p.Val <> p.val <> "+1"
EM_results <- data.frame(GO.ID=EM_pathways[,1], Description=EM_pathways[,1], 	
				p.Val=meta_top_pathways$Fisher_p, FDR=meta_top_pathways$Fisher_p, Phenotype="+1")		
EM_results[,3] <- format(EM_results[,3], scientific=FALSE)
EM_results[,4] <- format(EM_results[,4], scientific=FALSE)

# Write
setwd(output_dir)
write.table(EM_pathways, 'EM_pathways_fisher_untranslated.gmt', append=F, quote=F, row.names=F, col.names=F, na="", sep='\t')
write.table(EM_results, 'EM_results_fisher_untranslated.txt', append=F, quote=F, row.names=F, col.names=T, sep='\t')
##################################################################


##################################################################
##################################################################
# Put the two files above into Enrichment Map to make the network figure.
##################################################################
##################################################################



##################################################################
##################################################################
# Make heatmaps for GBJ
library(ggplot2)
library(viridis)
library(reshape2)

# Read merged data
setwd(output_dir)
results <- read.table('GBJ_merged_results.txt', header=T)
case_control_tab <- read.table('case_control_counts.txt', header=T)

# Add (n=) to the colnames of the results table
pvalue_tab <- results
rank_tab <- results
for (i in 3:ncol(pvalue_tab)) {
	
	# Build the rankings
	rank_tab[, i] <- rank(rank_tab[, i])
	
	# Add (n=)
	temp_study <- colnames(results[i])
	temp_row <- which(case_control_tab$Study == temp_study)
	temp_n <- case_control_tab$Cases[temp_row] +  case_control_tab$Controls[temp_row] 
	colnames(pvalue_tab)[i] <- paste(colnames(results[i]), ' (n=', temp_n, ')', sep='')
	colnames(rank_tab)[i] <- paste(colnames(results[i]), ' (n=', temp_n, ')', sep='')
}

# Melt both rankings and p-values
pvalue_melted <- melt(pvalue_tab, id.vars=c('Pathway', 'Num_genes'))
colnames(pvalue_melted) <- c('Pathway', 'Num_genes', 'Study', 'Pvalue')
pvalue_melted$Pvalue <- -log(pvalue_melted$Pvalue, base=10) 
rank_melted <- melt(rank_tab, id.vars=c('Pathway', 'Num_genes'))
colnames(rank_melted) <- c('Pathway', 'Num_genes', 'Study', 'Rank')

# Remove NAs
NA_rows <- which(is.na(pvalue_melted$Pvalue))
if( length(NA_rows) > 0) {
	pvalue_melted <- pvalue_melted[-NA_rows, ]
	rank_melted <- rank_melted[-NA_rows, ]
}

# Order y-axis of heatmap by sample size
pvalue_melted$Study <- factor(pvalue_melted$Study, levels=levels(pvalue_melted$Study)[order(rowSums(case_control_tab[, 2:3]))])
rank_melted$Study <- factor(rank_melted$Study, levels=levels(rank_melted$Study)[order(rowSums(case_control_tab[, 2:3]))])

# Order the x-axis of heatmap by the significance of METABRIC
ref_dat <- pvalue_melted[which(pvalue_melted$Study=='METABRIC (n=738)'), ]
ref_dat <- ref_dat[order(ref_dat$Pvalue, decreasing=TRUE),]
pvalue_melted$Pathway <- factor(pvalue_melted$Pathway, levels=ref_dat$Pathway)
rank_melted$Pathway <- factor(rank_melted$Pathway, levels=ref_dat$Pathway)

# Plot and save pvalue heatmap
#setEPS()
#postscript("SuppFigure3.eps")
jpeg('SuppFigure3.jpg')
legend_lab <- expression(paste(-log[10], "(p-value)", sep=''))
ggplot(pvalue_melted, aes(Pathway, Study, fill=Pvalue)) +
  geom_tile(color="white",size=0.1) +
  scale_fill_viridis(name=legend_lab) + 
  labs(x="Gene Sets (593 Total)", y="Dataset (21 Total)") +
  theme(axis.ticks=element_blank(), axis.text.x=element_blank()) 
dev.off()

# Plot and save rank heatmap
#setEPS()
#postscript("Figure2.eps")
jpeg('Figure2.jpg')
legend_lab <- expression(paste("Rank", sep=''))
ggplot(rank_melted, aes(Pathway, Study, fill=Rank)) +
  geom_tile(color="white",size=0.1) +
  scale_fill_viridis(name=legend_lab) + 
  labs(x="Gene Sets (593 Total)", y="Dataset (21 Total)") +
  theme(axis.ticks=element_blank(), axis.text.x=element_blank()) 
dev.off()




####################################################
# Unique SNP Sets
library(ggplot2)
setwd(ref_dir)

# Read merged data
GBJ_results <- read.table('GBJ_merged_results.txt', header=T)
GBJ_results <- GBJ_results[, -2]
GBJ_results <- GBJ_results[, -which(colnames(GBJ_results) %in% c('NCI', 'NKI', 'UCSF'))]
GSEA_results <- read.table('GSEA_merged_results.txt', header=T)
GSEA_results <- GSEA_results[, -2]
case_control_tab <- read.table('case_control_counts.txt', header=T)

# The number of different pathways found in the first x studies
threshold <- 59
num_unique_pathways <- data.frame(X=rep(1:threshold, 2), Test=rep(c('GBJ', 'GSEA'), each=threshold), Num_unique=NA)
for (i in 1:threshold) {
	
	# Clear list of names
	cat(i)
	unique_pathways_GBJ <- NULL
	unique_pathways_GSEA <- NULL
	
	# Go study by study, take the top x
	for (j in 2:ncol(GBJ_results)) {
		temp_GBJ_data <- GBJ_results[, c(1,j)]
		temp_GBJ_data <- temp_GBJ_data[order(temp_GBJ_data[,2]), ]
		temp_GBJ_data <- temp_GBJ_data[1:i, ]
		temp_GSEA_data <- GSEA_results[, c(1,j)]
		temp_GSEA_data <- temp_GSEA_data[order(temp_GSEA_data[,2]), ]
		temp_GSEA_data <- temp_GSEA_data[1:i, ]
		
		# Add to vector of names
		if (is.null(unique_pathways_GBJ)) {
			unique_pathways_GBJ <- temp_GBJ_data[, 1]
			unique_pathways_GSEA <- temp_GSEA_data[, 1]
		} else {
			unique_pathways_GBJ <- c(unique_pathways_GBJ, temp_GBJ_data[, 1])
			unique_pathways_GSEA <- c(unique_pathways_GSEA, temp_GSEA_data[, 1])
		}
	}
	
	# Record
	num_unique_pathways$Num_unique[i] <- length(unique(unique_pathways_GBJ))
	num_unique_pathways$Num_unique[i+threshold] <- length(unique(unique_pathways_GSEA))
}

# Plot
ggplot(data=num_unique_pathways, aes(x=X, y=Num_unique, group=Test)) + 
  geom_line(aes(color=Test)) + theme_bw() +
  scale_color_manual(values=c("darkorange", "darkblue")) +
  xlab('Number of sets reported from each study') + 
  ylab('Number of unique sets') + 
  theme(legend.position = c(0.9, 0.2), 
        axis.title=element_text(size=14, face='bold'),
        plot.title=element_text(size=14, face='bold'))
ggsave('Figure3.eps', device="eps")


####################################################
# Unique SNP Sets in small studies
library(ggplot2)
setwd(ref_dir)

# Read merged data
GBJ_results <- read.table('GBJ_merged_results.txt', header=T)
GBJ_results <- GBJ_results[, -2]
GBJ_results <- GBJ_results[, -which(colnames(GBJ_results) %in% c('NCI', 'NKI', 'UCSF'))]
GSEA_results <- read.table('GSEA_merged_results.txt', header=T)
GSEA_results <- GSEA_results[, -2]
case_control_tab <- read.table('case_control_counts.txt', header=T)

#Extract only 10 smallest studies
case_control_tab$Total <- case_control_tab$Cases + case_control_tab$Controls
gbj_studies <- which(names(GBJ_results) %in% case_control_tab[order(case_control_tab$Total),]$Study[1:12])
#Note this yields 10 studies
GBJ_results <- GBJ_results[,gbj_studies]
gsea_studies <- which(names(GSEA_results) %in% case_control_tab[order(case_control_tab$Total),]$Study[1:12])
GSEA_results <- GSEA_results[,gsea_studies]

# The number of different pathways found in the first x studies
threshold <- 59
num_unique_pathways <- data.frame(X=rep(1:threshold, 2), Test=rep(c('GBJ', 'GSEA'), each=threshold), Num_unique=NA)
for (i in 1:threshold) {
  
  # Clear list of names
  cat(i)
  unique_pathways_GBJ <- NULL
  unique_pathways_GSEA <- NULL
  
  # Go study by study, take the top x
  for (j in 2:ncol(GBJ_results)) {
    temp_GBJ_data <- GBJ_results[, c(1,j)]
    temp_GBJ_data <- temp_GBJ_data[order(temp_GBJ_data[,2]), ]
    temp_GBJ_data <- temp_GBJ_data[1:i, ]
    temp_GSEA_data <- GSEA_results[, c(1,j)]
    temp_GSEA_data <- temp_GSEA_data[order(temp_GSEA_data[,2]), ]
    temp_GSEA_data <- temp_GSEA_data[1:i, ]
    
    # Add to vector of names
    if (is.null(unique_pathways_GBJ)) {
      unique_pathways_GBJ <- temp_GBJ_data[, 1]
      unique_pathways_GSEA <- temp_GSEA_data[, 1]
    } else {
      unique_pathways_GBJ <- c(unique_pathways_GBJ, temp_GBJ_data[, 1])
      unique_pathways_GSEA <- c(unique_pathways_GSEA, temp_GSEA_data[, 1])
    }
  }
  
  # Record
  num_unique_pathways$Num_unique[i] <- length(unique(unique_pathways_GBJ))
  num_unique_pathways$Num_unique[i+threshold] <- length(unique(unique_pathways_GSEA))
}

# Plot
ggplot(data=num_unique_pathways, aes(x=X, y=Num_unique, group=Test)) + 
  geom_line(aes(color=Test)) + theme_bw() +
  scale_color_manual(values=c("darkorange", "darkblue")) +
  xlab('Number of sets reported from each study') + 
  ylab('Number of unique sets') + 
  theme(legend.position = c(0.9, 0.2), 
        axis.title=element_text(size=14, face='bold'),
        plot.title=element_text(size=14, face='bold'))
ggsave('SuppFigure4.eps', device="eps")




##################################################################
#Supplementary figure histograms
##################################################################
tfSets <- read.csv("tfSets.csv")
#Remove annotations
tfSetSubset <- tfSets[,3:dim(tfSets)[2]]
#Get count of number of genes in the set
geneCount <- apply(tfSetSubset, 1, function(x) sum(is.na(x)==FALSE))
library(ggplot2)
qplot(geneCount,
    geom="histogram",
    xlab = "Number of genes in TF set",
    ylab = "Count of TF sets",
    fill=I("deepskyblue4")) + theme_bw() +
    theme(axis.title=element_text(size=14, face='bold'))
ggsave('SuppFigure1.eps', device="eps")
median(geneCount)


#Get count of gene appearances (# sets a gene belongs to)
geneVec <- unlist(tfSetSubset)
geneVector <- as.character(geneVec[!is.na(geneVec)])
geneNums <- as.numeric(table(geneVector))
qplot(geneNums,
    geom="histogram",
    xlab = "Number of TF sets gene belongs to",
    ylab = "Count of genes ",
    fill=I("red4")) + theme_bw() +
    theme(axis.title=element_text(size=14, face='bold'))
ggsave('SuppFigure2.eps', device="eps")
median(table(geneVector))


	