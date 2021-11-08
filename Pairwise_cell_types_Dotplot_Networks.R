source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/Colour_Scheme.R")
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/SubColour_Scheme.R")
require(Seurat)



# Read in Data cell-type 1
express_dir <- "C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster"
type1 = "NKT"
obj1 <- readRDS(paste(express_dir,"/", type1, "_varimax_Subcluster.rds", sep = ""))
metadata1 <- readRDS(paste(express_dir,"/", type1, "_fullmetadata.rds", sep = ""))
subtypes1 <- unique(metadata1$Subcluster_Manual)


# Read in Data cell-type 2
express_dir <- "C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster"
type2 = "Endo"
obj2 <- readRDS(paste(express_dir,"/", type2, "_varimax_Subcluster.rds", sep = ""))
metadata2 <- readRDS(paste(express_dir,"/", type2, "_fullmetadata.rds", sep = ""))
subtypes2 <- unique(metadata2$Subcluster_Manual)

## Create Block list ##
all_genes <- unique(rownames(obj1), rownames(obj2))
gene_block_list <- c("SAA1", "ALB", all_genes[grepl("^RPS", all_genes)], all_genes[grepl("^RPL", all_genes)], all_genes[grepl("^MT-", all_genes)])

# Read in gene-gene interactions from cell-phonedb
files <- Sys.glob("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/CellPhonedb/*_cellphonedb_means.txt")
files <- files[!grepl("Flush", files)]

clean_up <- function(out, type1, type2) {
	# Clean up raw output from cellphonedb
	ids <- paste(out[,1], out[,2], sep="_")
	out[,1] <- ids;
	out <- out[,!grepl("Contam", colnames(out))]
	out <- out[,!grepl("Debris", colnames(out))]
	out <- out[,!grepl("Doublet", colnames(out))]
	out <- out[,c(1:11,which(grepl(type1, colnames(out)) & grepl(type2, colnames(out))))]
	return(out)
}

# Read in all output files
all_mats <- list()
all_comparisons <- c()
all_interactors <- c()
for (f in files) {
	mat <- read.table(f, header=T, sep="\t")
	mat <- clean_up(mat, type1, type2)
	tag <- sub("_cellphonedb_means.txt", "", f)
	all_mats[[tag]] <- mat
	
	all_interactors <- unique(rbind(all_interactors, mat[,1:11]))
	all_comparisons <- unique(c(all_comparisons, colnames(mat)[12:ncol(mat)]))
}
all_comparisons <- sort(all_comparisons)


# Loop over all sample matrices and add up the interaction scores for each pairwise comparison
	# selecting only ligand-receptor interactions.
align_matrix <- function(mat, rows, cols) {
	mat <- mat[ match(rows, mat[,1]) ,];
	tmp <- match(cols, colnames(mat))
	missing <- which(is.na(tmp));
	tmp[is.na(tmp)] <- 1

	mat <- mat[, tmp]
	mat[is.na(mat)] <- 0;
	mat[,missing] <- 0
	colnames(mat) <- cols;
	return(mat)
}

	
ligand_receptor_all <- all_interactors[ all_interactors[,"receptor_a"] != all_interactors[,"receptor_b"] & all_interactors[,"secreted"] == "False" ,]
total <- align_matrix(all_mats[[1]], ligand_receptor_all[,1], all_comparisons)
for( i in 2:length(all_mats)) {
	this_mat <- align_matrix(all_mats[[i]], ligand_receptor_all[,1], all_comparisons)
	total <- total+this_mat
}


# Then apply some threshold to ID the key genes for each pairwise comparison
threshold <- signif(quantile( apply(total, 1, max), 1-40/nrow(total)), digits=1)
#threshold <- 5  # 7 = 0.5*length(all_mats)/2

# Then collect these genes sorting them by ligand or receprot
keep <- rowSums(total > threshold) > 0
ligands <- unique(c(ligand_receptor_all[keep & ligand_receptor_all[,"receptor_a"] == "False", "gene_a"],
		ligand_receptor_all[keep & ligand_receptor_all[,"receptor_b"] == "False", "gene_b"]))
receptors <- unique(c(ligand_receptor_all[keep & ligand_receptor_all[,"receptor_a"] == "True", "gene_a"],
		ligand_receptor_all[keep & ligand_receptor_all[,"receptor_b"] == "True", "gene_b"]))

# Collect top ligand-receptor pairs for this pair of cell-types
genes1 <- c(ligands, receptors)
genes2 <- c(receptors, ligands)
genes1 <- genes1[!genes1 %in% gene_block_list]
genes2 <- genes2[!genes2 %in% gene_block_list]

# Plot Dotplot 1
require(ggplot2)
obj1@meta.data <- metadata1
obj1_clean <- obj1[,!grepl("Contam", obj1@meta.data$Subcluster_Manual)]
obj1_clean <- obj1_clean[,!grepl("Debris", obj1_clean@meta.data$Subcluster_Manual)]
obj1_clean <- obj1_clean[,!grepl("Doublet", obj1_clean@meta.data$Subcluster_Manual)]
obj1_clean <- obj1_clean[,!grepl("Flush", obj1_clean@meta.data$Subcluster_Manual)]
obj1_clean <- obj1_clean[,!grepl("Hepato", obj1_clean@meta.data$Subcluster_Manual)]

#Detection Rate - Keep only genes detected in at least 10% of cells of one of the clusters
my_group_rowmeans <- function(MAT, group_labs, type=c("mean","sum")) {
        d <- split(1:ncol(MAT), group_labs);
        if (type[1] == "mean") {
                mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
        } else {
                mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
        } 
        return(mus);
}

detect <- my_group_rowmeans(obj1_clean@assays$RNA@counts[rownames(obj1_clean) %in% genes1,] > 0, factor(obj1_clean@meta.data$Subcluster_Manual))
detect <- detect[match(genes1, rownames(detect)),]
genes1 <- genes1[apply(detect, 1, max) > 0.15 & !is.na(rowSums(detect))]

DotPlot(obj1_clean, features=genes1, group.by="Subcluster_Manual")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot Dotplot 2
require(ggplot2)

obj2@meta.data <- metadata2
obj2_clean <- obj2[,!grepl("Contam", obj2@meta.data$Subcluster_Manual)]
obj2_clean <- obj2_clean[,!grepl("Debris", obj2_clean@meta.data$Subcluster_Manual)]
obj2_clean <- obj2_clean[,!grepl("Doublet", obj2_clean@meta.data$Subcluster_Manual)]
obj2_clean <- obj2_clean[,!grepl("Flush", obj2_clean@meta.data$Subcluster_Manual)]
obj2_clean <- obj2_clean[,!grepl("Hepato", obj2_clean@meta.data$Subcluster_Manual)]


#Detection Rate - Keep only genes detected in at least 10% of cells of one of the clusters
my_group_rowmeans <- function(MAT, group_labs, type=c("mean","sum")) {
        d <- split(1:ncol(MAT), group_labs);
        if (type[1] == "mean") {
                mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
        } else {
                mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
        } 
        return(mus);
}

detect <- my_group_rowmeans(obj2_clean@assays$RNA@counts[rownames(obj2_clean) %in% genes2,] > 0, factor(obj2_clean@meta.data$Subcluster_Manual))
detect <- detect[match(genes2, rownames(detect)),]
genes2 <- genes2[apply(detect, 1, max) > 0.15 & !is.na(rowSums(detect))]

DotPlot(obj2_clean, features=genes2, group.by="Subcluster_Manual")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot connections between them
arrows_to_make <- ligand_receptor_all[keep,]
arrows_from <- c(arrows_to_make[ arrows_to_make[,"receptor_a"] == "False" ,"gene_a"], 
			arrows_to_make[ arrows_to_make[,"receptor_b"] == "False" ,"gene_b"])

arrows_to <- c(arrows_to_make[ arrows_to_make[,"receptor_a"] == "True" ,"gene_a"], 
			arrows_to_make[ arrows_to_make[,"receptor_b"] == "True" ,"gene_b"])

cbind(arrows_from, arrows_to)
