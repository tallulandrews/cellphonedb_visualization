# Combine all mean files into one.

# Notes: Cellphonedb output is directional but this directionality is arbitrary : simply 
# based on which protein is labeled A or B for each interating pair. Therefore, I shoul
# make the cell-cell interaction matrix undirectional.

# Could also subset to just interactions that have directionality and rearrange matrix accordingly.

dir <- "/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Subcluster/ManualAnnotation"


files <- Sys.glob("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/CellPhonedb/*_cellphonedb_means.txt")
files <- files[!grepl("Flush", files)]

clean_up <- function(out) {
	# Clean up raw output from cellphonedb
	ids <- paste(out[,1], out[,2], sep="_")
	out[,1] <- ids;
	out <- out[,!grepl("Contam", colnames(out))]
	out <- out[,!grepl("Debris", colnames(out))]
	out <- out[,!grepl("Doublet", colnames(out))]
	return(out)
}

# Read in all output files
all_mats <- list()
all_comparisons <- c()
all_interactors <- c()
for (f in files) {
	mat <- read.table(f, header=T, sep="\t")
	mat <- clean_up(mat)
	tag <- sub("_cellphonedb_means.txt", "", f)
	all_mats[[tag]] <- mat
	
	all_interactors <- unique(rbind(all_interactors, mat[,1:11]))
	all_comparisons <- unique(c(all_comparisons, colnames(mat)[12:ncol(mat)]))
}

# Thresholds for what is considered a valid interaction:

expr_threshold = 0.5; # expression level required of the interaction to keep it
n_interaction_threshold = 5; # number of interacting pairs to declare two cell-type interacting.

all_comparisons <- sort(all_comparisons)


# Count number ofsample in which each pair of cell-types are "interacting" 
# threshold of expression for a pair of proteins to be interacting
# threshold of number of pairs for a cell-type to be interating
significance_n_edges <- function(mat, expr_threshold=0.5) {
	mat <- ( mat+t(mat) )/2
	n_interpairs <- colSums(mat > expr_threshold, na.rm=T)
	

}

# Given a list of matrices for integer values of how many interactions are between A & B
# Where each matrix is from a different sample.
# Is there more interaction across multiple indivuals than we expect by chance?
# Conrolling for the fact that certain cell-types express more receptors/ligands than others.

# Chi Square test?(O-E)^2/E
# Expected number of interactions between A & B = sum_A(edges)/nodes * sum_B(edges)/nodes
# Dot plot? Average obs-expect & number of samples where > 0?
Obs_vs_Expect_nEdges <- function(list_of_mats) {
	
	mat <- list_of_mats[[i]]
	expected <-  (rowSums(mat) %*% t(colSums(mat)))/sum(mat);
	mat-expected

}

# if expr_threshold < 0 uses sum of expression rather than number of interactions.
convert_interactionMat_2_type_by_type <- function(out, expr_threshold=0.5) {
	types <- get_cell_types(colnames(out))
	all_types <- unique(types[,1])	
	if (expr_threshold < 0) {
		n_interaction <- colSums(out, na.rm=T)
	} else {
		is.interaction <- out > expr_threshold
		n_interaction <- colSums(is.interaction, na.rm=T)
	}

	type_by_type <- matrix(0, nrow=length(all_types), ncol=length(all_types));
	for(r in 1:nrow(types)) {
		i <- which((all_types) == types[r,1])
		j <- which((all_types) == types[r,2])
		type_by_type[i,j] <- n_interaction[r]
	}
	colnames(type_by_type) <- all_types
	rownames(type_by_type) <- all_types
	return(type_by_type)
}

# Parse cellphonedb column labels into cell-type pairs.
get_cell_types <- function(comparisons) {
	comparisons <- gsub("\\.qHSC", "_qHSC", comparisons)
	comparisons <- gsub("\\.aHSC", "_aHSC", comparisons)
	cell_types <- strsplit(comparisons, "\\.")
	type_a <- sapply(cell_types, function(x){unlist(x)[1]})
	type_b <- sapply(cell_types, function(x){unlist(x)[min(length(x),2)]})
	return(cbind(type_a, type_b))
}
	

convert_cellphonedb_2_ligand_receptor <- function(out, expr_threshold=0.5) {
	# Receptor = a = first type = row
	# Ligand = b = second type = column
	out <- out[ out[, "secreted"]=="False" & out[, "receptor_a"] != out[, "receptor_b"],]

	cell_types <- get_cell_types(colnames(out))
	type_a <- cell_types[,1]
	type_b <- cell_types[,2]
	for (i in which(out[,"receptor_a"] == "False")) {
		for (j in ( which(out[i,12:ncol(out)] > 0)+11 ) ){
			if (type_a[j] == type_b[j]) {
				next;
			}
			swap <- which(type_b == type_a[j] & type_a == type_b[j])
			tmp <- out[i,swap]
			out[i,swap] <- out[i,j]
			out[i,j] <- tmp
		}
	}
	mat <- convert_interactionMat_2_type_by_type(out[,12:ncol(out)], expr_threshold=expr_threshold)
	return(mat)
}

convert_cellphonedb_2_secreted_receptor <- function(out, expr_threshold=0.5) {
	# Receptor = a = first type = row
	# Ligand = b = second type = column
	out <- out[ out[, "secreted"]=="True" & out[, "receptor_a"] != out[, "receptor_b"],]

	cell_types <- get_cell_types(colnames(out))
	type_a <- cell_types[,1]
	type_b <- cell_types[,2]
	for (i in which(out[,"receptor_a"] == "False")) {
		for (j in ( which(out[i,12:ncol(out)] > 0)+11 ) ){
			if (type_a[j] == type_b[j]) {
				next;
			}
			swap <- which(type_b == type_a[j] & type_a == type_b[j])
			tmp <- out[i,swap]
			out[i,swap] <- out[i,j]
			out[i,j] <- tmp
		}
	}
	mat <- convert_interactionMat_2_type_by_type(out[,12:ncol(out)], expr_threshold=expr_threshold)
	return(mat)
}


convert_cellphonedbe_2_local <- function(out, expr_threshold=0) {
        # Receptor = a = first type = row
        # Ligand = b = second type = column
        out <- out[ out[, "secreted"]=="False",]

        mat <- convert_interactionMat_2_type_by_type(out[,12:ncol(out)], expr_threshold=expr_threshold)
	mat <- mat + t(mat)
	return(mat)
}


###### Making Plots #######
names(all_mats) <- lapply(strsplit(names(all_mats),"/"), function(x) {x[[length(x)]]})


### Preliminary Heatmaps  - overall strength of any kind of interactions ###
full_matrix <- c()
n_mat <- c();
sum_matrix <- c();
full_matrix_cell_types <- c()
full_matrix_sample <- c();



n_present <- rep(0, length(all_comparisons)) # 
n_sample_inter <- rep(0, length(all_comparisons)) # Number of samples where a particular pair are interacting

# Each sample separately

for (i in 1:length(all_mats)) {
	this_mat <- all_mats[[i]]
	this_mat <- this_mat[match( all_interactors[,1], this_mat[,1]),] # Fix cases of missing cell-types in some samples, other missing genes
	rownames(this_mat) <- all_interactors[,1]
	tmp <- as.matrix(this_mat[,colnames(this_mat) %in% all_comparisons])
	tmp <- tmp[,match(all_comparisons, colnames(tmp))]
	colnames(tmp) <- paste(names(all_mats)[i], all_comparisons, sep="_")
	this_mat <- tmp;
	full_matrix_cell_types <- c(full_matrix_cell_types, all_comparisons)
	full_matrix_sample <- c(full_matrix_sample, rep(names(all_mats)[i], length(all_comparisons)))

	if (is.null(nrow(full_matrix))) {
		full_matrix <- this_mat
		this_mat[is.na(this_mat)] <- 0;
		sum_matrix <- this_mat;
		n_mat <- !is.na(this_mat) + 0;
	} else {
		full_matrix <- cbind(full_matrix, this_mat)
		n_mat <- n_mat + !is.na(this_mat);
		this_mat[is.na(this_mat)] <- 0;
		sum_matrix <- sum_matrix + this_mat;
	}

	# Count T/F interactions
	present <- !( colSums(is.na(tmp)) == nrow(tmp))
	n_present <- n_present + present
	n_interpairs <- colSums(tmp > expr_threshold, na.rm=T)
	n_sample_inter <- n_sample_inter + (n_interpairs > n_interaction_threshold)
}

mean_across_samples_mat <- sum_matrix/n_mat
colnames(mean_across_samples_mat) <- lapply(strsplit(colnames(mean_across_samples_mat),"/"), function(x) {x[[length(x)]]})

overall_score <- colSums(mean_across_samples_mat)
cell_types <- sub("C37_", "", names(overall_score))
cell_types <- gsub("\\.qHSC", "_qHSC", cell_types)
cell_types <- gsub("\\.aHSC", "_aHSC", cell_types)

cell_types <- strsplit(cell_types, "\\.")

cell_type_2_id <- factor(unique(unlist(cell_types))) ### Very Important for Later!!!

heatmap_dat <- matrix(0, nrow=length(cell_type_2_id), ncol=length(cell_type_2_id))
heatmap_nsamp_dat <- matrix(0, nrow=length(cell_type_2_id), ncol=length(cell_type_2_id))
for(r in 1:length(cell_types)) {
	i <- which(levels(cell_type_2_id) == cell_types[[r]][1])
	j <- which(levels(cell_type_2_id) == cell_types[[r]][2])
	heatmap_dat[i,j] <- overall_score[r]
	heatmap_nsamp_dat[i,j] <- n_sample_inter[r]/n_present[r]
}

heatmap_dat <- (heatmap_dat + t(heatmap_dat))/2
colnames(heatmap_dat) <- rownames(heatmap_dat) <- levels(cell_type_2_id)
heatmap_nsamp_dat <- (heatmap_nsamp_dat + t(heatmap_nsamp_dat))/2
colnames(heatmap_nsamp_dat) <- rownames(heatmap_nsamp_dat) <- levels(cell_type_2_id)

require(gplots)

png("Mean_pair_expr_prelim.png", width=10, height=10, units="in", res=300)
heatmap.2(heatmap_dat, margins=c(8, 8), trace="none")
dev.off()
png("Mean_pair_replicated_prelim.png", width=10, height=10, units="in", res=300)
heatmap.2(heatmap_nsamp_dat, margins=c(8, 8), trace="none")
dev.off()

# Observed vs Expected #
require(RColorBrewer)
require(Seurat)
palette <- PurpleAndYellow(20)

png("Mean_pair_expr_vs_expected_prelim.png", width=10, height=10, units="in", res=300)
expected <-  (rowSums(heatmap_dat) %*% t(colSums(heatmap_dat)))/sum(heatmap_dat);
heatmap.2(heatmap_dat-expected, margins=c(8, 8), trace="none", col=palette)
dev.off()

png("Mean_pair_replicated_vs_expected_prelim.png", width=10, height=10, units="in", res=300)
expected <-  (rowSums(heatmap_nsamp_dat) %*% t(colSums(heatmap_nsamp_dat)))/sum(heatmap_nsamp_dat);
heatmap.2(heatmap_nsamp_dat-expected, margins=c(8, 8), trace="none", col=palette)
dev.off()

######## Ligand-Receptor - Surface Bound #####
# Requires from Above: cell_type_2_id
n_mat <- matrix(0, ncol=length(cell_type_2_id), nrow=length(cell_type_2_id)); # Number of samples where samples both cell_types are present so could assay their interactions.
n_interact <- matrix(0, ncol=length(cell_type_2_id), nrow=length(cell_type_2_id)); # Number of samples where >= 1 interacting pair.
sum_matrix <- matrix(0, ncol=length(cell_type_2_id), nrow=length(cell_type_2_id)); # Total Interactions across all samples.

rownames(sum_matrix) <- cell_type_2_id; colnames(sum_matrix) <- cell_type_2_id;
rownames(n_mat) <- cell_type_2_id; colnames(n_mat) <- cell_type_2_id;
rownames(n_interact) <- cell_type_2_id; colnames(n_interact) <- cell_type_2_id;


# Each sample separately

for (i in 1:length(all_mats)) {
	this_mat <- all_mats[[i]]
	this_ligand_receptor <- convert_cellphonedb_2_ligand_receptor(this_mat, expr_threshold=0.5)

	this_ligand_receptor <- this_ligand_receptor[match(as.character(cell_type_2_id), rownames(this_ligand_receptor)),
									match(as.character(cell_type_2_id), colnames(this_ligand_receptor))]

	n_mat <- n_mat + !is.na(this_ligand_receptor);
	this_ligand_receptor[is.na(this_ligand_receptor)] <- 0
	sum_matrix <- sum_matrix+this_ligand_receptor
	n_interact <- n_interact+( this_ligand_receptor > 0 )
}


require(gplots)


png("Receptor_by_Ligand_mean.png", width=10, height=10, units="in", res=300)
heatmap_dat <- sum_matrix/n_mat;
heatmap.2(heatmap_dat, margins=c(8, 8), trace="none")
dev.off()


## For Job Talk ##

colours <- c("grey65", "#fdbf6f", "#ff7f00", "#ffff99", "#a6cee3", "#1f78b4")
exclude <- grep("Hepatocyte2", colnames(heatmap_dat))
this_dat <- heatmap_dat[-exclude, -exclude]
exclude <- which(colnames(this_dat) %in% c("NKT_Erythroblasts", "NKT_Prolif",
		"NKT_Hepato", "ProlifInterzonalHep", "NKT_Flush", "Stellate_aHSC",
		"Stellate_AP1_", "Stellate_Fibro_aHSC", "NKT_Myeloid", "Macrophage_RetNonInf"))
this_dat <- this_dat[-exclude, -exclude]
colnames(this_dat) <- rownames(this_dat) <- c("Bcell_IgKA", "Bcell_IgKG", "Bcell_IgLA", "Bcell_IgLG",
		"Bcell_Naive", "Bcell_Prolif", "HepCentral", "Chol", "Chol2", "Chol_Mucus",
		"cvEndo", "cvLSEC", "interLSEC", "ppLSEC", "VasEndo", "HepInter",
		"aInf_Mac", "Inf_Mac", "MHCII_Mac", "NonInf_Mac", "Phago_Mac", 
		"aNonInf_Mac", "Synap_Mac", "CD3_Tcell", "CD8_Tcell", "cNKcell", "lrNKcell", "Bcell_Mature",
		"HepPortal", "Stellate_Fibro", "Stellate_qHSC")



require(Seurat)
png("Receptor_by_Ligand_mean_cleanedup.png", width=10, height=10, units="in", res=300)
heatmap.2(this_dat, col=CustomPalette(low="white", high="black", k=20), trace="none", margins=c(8, 8))
dev.off()

# Observed vs Expected #
require(RColorBrewer)
require(Seurat)
palette <- PurpleAndYellow(20)

png("Receptor_by_Ligand_mean_vs_expected.png", width=10, height=10, units="in", res=300)
expected <-  (rowSums(heatmap_dat) %*% t(colSums(heatmap_dat)))/sum(heatmap_dat);
heatmap.2(heatmap_dat-expected, margins=c(8, 8), trace="none", col=palette)
dev.off()

######## Ligand-Receptor - Secreted #####
# Requires from Above: cell_type_2_id
n_mat <- matrix(0, ncol=length(cell_type_2_id), nrow=length(cell_type_2_id)); # Number of samples where samples both cell_types are present so could assay their interactions.
n_interact <- matrix(0, ncol=length(cell_type_2_id), nrow=length(cell_type_2_id)); # Number of samples where >= 1 interacting pair.
sum_matrix <- matrix(0, ncol=length(cell_type_2_id), nrow=length(cell_type_2_id)); # Total Interactions across all samples.

rownames(sum_matrix) <- cell_type_2_id; colnames(sum_matrix) <- cell_type_2_id;
rownames(n_mat) <- cell_type_2_id; colnames(n_mat) <- cell_type_2_id;
rownames(n_interact) <- cell_type_2_id; colnames(n_interact) <- cell_type_2_id;

# Each sample separately

for (i in 1:length(all_mats)) {
	this_mat <- all_mats[[i]]
	this_ligand_receptor <- convert_cellphonedb_2_secreted_receptor(this_mat, expr_threshold=0.5)

	this_ligand_receptor <- this_ligand_receptor[match(as.character(cell_type_2_id), rownames(this_ligand_receptor)),
									match(as.character(cell_type_2_id), colnames(this_ligand_receptor))]

	n_mat <- n_mat + !is.na(this_ligand_receptor);
	this_ligand_receptor[is.na(this_ligand_receptor)] <- 0
	sum_matrix <- sum_matrix+this_ligand_receptor
	n_interact <- n_interact+( this_ligand_receptor > 0 )
}


require(gplots)

png("Receptor_by_Secreted_mean.png", width=10, height=10, units="in", res=300)
heatmap_dat <- sum_matrix/n_mat;
heatmap.2(heatmap_dat, margins=c(8, 8), trace="none")
dev.off()


# Observed vs Expected #
require(RColorBrewer)
require(Seurat)
palette <- PurpleAndYellow(20)

png("Receptor_by_Secreted_mean_vs_expected.png", width=10, height=10, units="in", res=300)
expected <-  (rowSums(heatmap_dat) %*% t(colSums(heatmap_dat)))/sum(heatmap_dat);
heatmap.2(heatmap_dat-expected, margins=c(8, 8), trace="none", col=palette)
dev.off()


## Macrophage Endothelial Ligand-Receptor ##

ligand_receptor <- all_interactors$receptor_a != all_interactors$receptor_b & all_interactors$secreted == "False"
Macro_Endo_columns <- grepl("Macrophage", colnames(full_matrix)) & grepl("Endo", colnames(full_matrix))

macro_endo <- full_matrix[ligand_receptor, Macro_Endo_columns]
score2 <- rowSums(macro_endo > 0, na.rm=TRUE)
keep <- score2 > 100
macro_endo <- macro_endo[keep,]
gene_info <- all_interactors[match(rownames(macro_endo), all_interactors[,1]),]

# Create edges / directionality
# Average strength of the interactions where Macrophage is A vs where Macrophage is B (across all pairs of interactions)
cell_types <- get_cell_types(colnames(macro_endo))
score_macro_to_endo <- rowMeans(macro_endo[, grepl("Macrophage", cell_types[,1]) & grepl("Endo", cell_types[,2])], na.rm=T)
score_endo_to_macro <- rowMeans(macro_endo[, grepl("Endo", cell_types[,1]) & grepl("Macro", cell_types[,2])], na.rm=T)

# Stronger interactions when one of the macrophage groups is A or when they are B?
is_a_macro <- score_macro_to_endo > score_endo_to_macro

# Get genes where Endo is ligand and Macro is receptor
A_is_macro <- gene_info[is_a_macro,]
B_is_macro <- gene_info[!is_a_macro,]
endo_to_macro <- rbind(A_is_macro[A_is_macro[,"receptor_a"] == "True",], B_is_macro[B_is_macro[,"receptor_b"] == "True",])

endo_genes <- c(A_is_macro[A_is_macro[,"receptor_a"] == "True","gene_b"], B_is_macro[B_is_macro[,"receptor_b"] == "True","gene_a"])
macro_genes <- c(A_is_macro[A_is_macro[,"receptor_a"] == "True","gene_a"], B_is_macro[B_is_macro[,"receptor_b"] == "True","gene_b"])

# Get expression of these genes across samples
macro_expression <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Macrophage_harmony_subcluster_pseudobulks.rds")
endo_expression <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Endo_harmony_subcluster_pseudobulks.rds")

macro_expr <- rowMeans(macro_expression$pseudobulk_mat[match(macro_genes,rownames(macro_expression$pseudobulk_mat)),])
endo_expr <- rowMeans(endo_expression$pseudobulk_mat[match(endo_genes,rownames(endo_expression$pseudobulk_mat)),])

keep <- ! (is.na(macro_expr) | is.na(endo_expr))
macro_expr <- macro_expr[keep]
endo_expr <- endo_expr[keep]

#Turn into graph
require(igraph)
G<-graph_from_edgelist(cbind(names(endo_expr), names(macro_expr)))

node_colour <- rep("#fdbf6f", length(V(G)))
node_colour[names(V(G)) %in% macro_genes] <- "#a6cee3"

png("Endothelial_to_Macrophage.png", width=7.5, height=9, units="in", res=300)
plot(G, edge.color="black", vertex.color=node_colour, label.cex=0.7, label.color="black", vertex.frame.color=NA, edge.arrow.size=0.5, edge.width=2)
dev.off()

## NKT Endothelial Ligand-Receptor ##

ligand_receptor <- all_interactors$receptor_a != all_interactors$receptor_b & all_interactors$secreted == "False"
NKT_Endo_columns <- grepl("NKT", colnames(full_matrix)) & grepl("Endo", colnames(full_matrix))

nkt_endo <- full_matrix[ligand_receptor, NKT_Endo_columns]
score2 <- rowSums(nkt_endo > 0, na.rm=TRUE)
keep <- score2 > 100
nkt_endo <- nkt_endo[keep,]
gene_info <- all_interactors[match(rownames(nkt_endo), all_interactors[,1]),]

# Create edges / directionality
# Average strength of the interactions where Macrophage is A vs where Macrophage is B (across all pairs of interactions)
cell_types <- get_cell_types(colnames(nkt_endo))
score_nkt_to_endo <- rowMeans(nkt_endo[, grepl("NKT", cell_types[,1]) & grepl("Endo", cell_types[,2])], na.rm=T)
score_endo_to_nkt <- rowMeans(nkt_endo[, grepl("Endo", cell_types[,1]) & grepl("NKT", cell_types[,2])], na.rm=T)

# Stronger interactions when one of the macrophage groups is A or when they are B?
is_a_nkt <- score_nkt_to_endo > score_endo_to_nkt

# Get genes where Endo is ligand and Macro is receptor
A_is_nkt <- gene_info[is_a_nkt,]
B_is_nkt <- gene_info[!is_a_nkt,]
endo_to_nkt <- rbind(A_is_nkt[A_is_nkt[,"receptor_a"] == "True",], B_is_nkt[B_is_nkt[,"receptor_b"] == "True",])

endo_genes <- c(A_is_nkt[A_is_nkt[,"receptor_a"] == "True","gene_b"], B_is_nkt[B_is_nkt[,"receptor_b"] == "True","gene_a"])
nkt_genes <- c(A_is_nkt[A_is_nkt[,"receptor_a"] == "True","gene_a"], B_is_nkt[B_is_nkt[,"receptor_b"] == "True","gene_b"])

# Get expression of these genes across samples
nkt_expression <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/NKT_harmony_subcluster_pseudobulks.rds")
endo_expression <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Endo_harmony_subcluster_pseudobulks.rds")

nkt_expr <- rowMeans(nkt_expression$pseudobulk_mat[match(nkt_genes,rownames(nkt_expression$pseudobulk_mat)),])
endo_expr <- rowMeans(endo_expression$pseudobulk_mat[match(endo_genes,rownames(endo_expression$pseudobulk_mat)),])

keep <- ! (is.na(nkt_expr) | is.na(endo_expr))
nkt_expr <- nkt_expr[keep]
endo_expr <- endo_expr[keep]

#Turn into graph
require(igraph)
G<-graph_from_edgelist(cbind(names(endo_expr), names(nkt_expr)))

node_colour <- rep("#fdbf6f", length(V(G)))
node_colour[names(V(G)) %in% nkt_genes] <- "#4daf4a"

png("Endothelial_to_NKT.png", width=7.5, height=9, units="in", res=300)
plot(G, edge.color="black", vertex.color=node_colour, label.cex=0.7, label.color="black", vertex.frame.color=NA, edge.arrow.size=0.5, edge.width=2)
dev.off()


###### Generalized Pairwise Cell-type interaction Plotting ######
interaction_filter <- all_interactors$receptor_a != all_interactors$receptor_b & all_interactors$secreted == "False"
cell_type_filter <- grepl("Macrophage", colnames(full_matrix)) & grepl("NKT", colnames(full_matrix))

# Colours: NKT = "#4daf4a", Mac = "#a6cee3", LSEC = "#fdbf6f"
cell_type_color_1 = "#a6cee3"
cell_type_color_2 = "#4daf4a"

## Gene expression in both cell-types - detection rate in each cell-type
type1_expression <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/Macrophage_harmony_Subcluster_Allgenes.rds")
type1_metadata <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Macrophage_fullmetadata.rds")
detect_3pr_1 <- group_rowmeans(type1_expression@assays$RNA@counts[, type1_expression@meta.data$assay_type == "3pr"] > 0, type1_metadata$Subcluster_Manual[type1_expression@meta.data$assay_type == "3pr"])
detect_5pr_1 <- group_rowmeans(type1_expression@assays$RNA@counts[, type1_expression@meta.data$assay_type == "5pr"] > 0, type1_metadata$Subcluster_Manual[type1_expression@meta.data$assay_type == "5pr"])
rm(type1_expression)

clean_detect_mat <- function(mat) {
	types <- colnames(mat)
	exclude <- grepl("Hepato", types) | grepl("Debris", types) | grepl("Doublet", types)
	mat <- mat[,!exclude]
	return(mat)
}
detect_3pr_1 <- rowMeans(clean_detect_mat(detect_3pr_1))
detect_5pr_1 <- rowMeans(clean_detect_mat(detect_5pr_1))


type2_expression <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/NKT_harmony_Subcluster_Allgenes.rds")
type2_metadata <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/NKT_fullmetadata.rds")
detect_3pr_2 <- group_rowmeans(type2_expression@assays$RNA@counts[, type2_expression@meta.data$assay_type == "3pr"] > 0, type2_metadata$Subcluster_Manual[type2_expression@meta.data$assay_type == "3pr"])
detect_5pr_2 <- group_rowmeans(type2_expression@assays$RNA@counts[, type2_expression@meta.data$assay_type == "5pr"] > 0, type2_metadata$Subcluster_Manual[type2_expression@meta.data$assay_type == "5pr"])
rm(type2_expression)

detect_3pr_2 <- rowMeans(clean_detect_mat(detect_3pr_2))
detect_5pr_2 <- rowMeans(clean_detect_mat(detect_5pr_2))
##

## Filter interactions
these_interactions <- full_matrix[interaction_filter, cell_type_filter]
score <- rowSums(these_interactions > 0, na.rm=TRUE)
keep <- score > rev(sort(score))[min(51, length(score))]; # top 50
these_interactions <- these_interactions[keep,]
gene_info <- all_interactors[match(rownames(these_interactions), all_interactors[,1]),]

## complexes to genes ##
complex_table <- read.table("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/CellPhonedb/CellPhonedb_database/complex_input.csv", sep=",", header=T)
gene_id_table <- read.table("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/CellPhonedb/CellPhonedb_database/gene_input.csv", sep=",", header=T)
map_id <- function(uniprot) {
	hgnc <- gene_id_table[match(uniprot, gene_id_table[,2]), 3]
	return(hgnc[!is.na(hgnc)])
}
complex_2_genes <- apply(complex_table, 1, function(x) { c( map_id( x[2:5] ) )})
names(complex_2_genes) <- complex_table[,1]
					
# Figure out directionality
fix_complexes <- function(genes) {
		fixed <- lapply(genes, function(x) {if (grepl("complex:",x)) {return(complex_2_genes[[sub("complex:", "", x)]])} else {return(x)}})
		names(fixed) <- genes
		return(fixed)
}

gene_1 <- gene_info[,5]; gene_1[gene_1 == ""] <- gene_info[gene_1 == "",3]; gene_1 <- fix_complexes(gene_1)
gene_2 <- gene_info[,6]; gene_2[gene_2 == ""] <- gene_info[gene_2 == "",4]; gene_2 <- fix_complexes(gene_2)

gene_1_3pr_diff <- sapply(gene_1, function(x) { mean(detect_3pr_1[names(detect_3pr_1) %in% x]) }) - sapply(gene_1, function(x) { mean(detect_3pr_2[names(detect_3pr_2) %in% x]) }) # > 0 = type 2 > type 1
gene_1_5pr_diff <- sapply(gene_1, function(x) { mean(detect_5pr_1[names(detect_5pr_1) %in% x]) }) - sapply(gene_1, function(x) { mean(detect_5pr_2[names(detect_5pr_2) %in% x]) })

gene_2_3pr_diff <- sapply(gene_2, function(x) { mean(detect_3pr_1[names(detect_3pr_1) %in% x]) }) - sapply(gene_2, function(x) { mean(detect_3pr_2[names(detect_3pr_2) %in% x]) })
gene_2_5pr_diff <- sapply(gene_2, function(x) { mean(detect_5pr_1[names(detect_5pr_1) %in% x]) }) - sapply(gene_2, function(x) { mean(detect_5pr_2[names(detect_5pr_2) %in% x]) })

expr_direction <- (gene_1_3pr_diff+gene_1_5pr_diff) - (gene_2_3pr_diff+gene_2_5pr_diff)
# expr_dir > 0 -> gene1: type2 > type1 or gene2: type1 > type2
# expr_dir < 0 -> gene1: type1 > type2 or gene2: type2 > type1
gene_1_colour <- rep(

interaction_direction <- rep(0, nrow(gene_info)) # 0 = both receptor / ligand, 1 = gene 1 is receptor, -1 = gene 2 is receptor
interaction_direction[gene_info[,8] == "True" & gene_info[,9] == "False"] = 1
interaction_direction[gene_info[,8] == "False" & gene_info[,9] == "True"] = -1

edge_direction <- sign(expr_direction)*interaction_direction
# edge_direction == 1 : type1 -> type2
# edge_direction == -1: type2 -> type1

# Nodes





#------------------------------------------------------------
# Average strength of the interactions where Macrophage is A vs where Macrophage is B (across all pairs of interactions)
cell_types <- get_cell_types(colnames(these_interaction))
score_nkt_to_endo <- rowMeans(nkt_endo[, grepl("NKT", cell_types[,1]) & grepl("Endo", cell_types[,2])], na.rm=T)
score_endo_to_nkt <- rowMeans(nkt_endo[, grepl("Endo", cell_types[,1]) & grepl("NKT", cell_types[,2])], na.rm=T)

# Stronger interactions when one of the macrophage groups is A or when they are B?
is_a_nkt <- score_nkt_to_endo > score_endo_to_nkt

# Get genes where Endo is ligand and Macro is receptor
A_is_nkt <- gene_info[is_a_nkt,]
B_is_nkt <- gene_info[!is_a_nkt,]
endo_to_nkt <- rbind(A_is_nkt[A_is_nkt[,"receptor_a"] == "True",], B_is_nkt[B_is_nkt[,"receptor_b"] == "True",])

endo_genes <- c(A_is_nkt[A_is_nkt[,"receptor_a"] == "True","gene_b"], B_is_nkt[B_is_nkt[,"receptor_b"] == "True","gene_a"])
nkt_genes <- c(A_is_nkt[A_is_nkt[,"receptor_a"] == "True","gene_a"], B_is_nkt[B_is_nkt[,"receptor_b"] == "True","gene_b"])





############ Graphs are unreadable ########################

## Global Graph plots
# Edge thickness = Number of interacting pairs (avg across samples)
# Edge Colour = expression level of interacting pairs (average across samples & pairs)
sum_exprLevel <- c();
sum_interPairs  <- c();
sum_present <- c(); # get this from n_presnt


# Graph 1 = Reciprocal any-any - Graphs are too messy! Stick to Heatmaps!

for (i in 1:length(all_mats)) {
        this_mat <- all_mats[[i]]
        this_mat <- this_mat[match( all_interactors[,1], this_mat[,1]),]
        rownames(this_mat) <- all_interactors[,1]
        tmp <- as.matrix(this_mat[,colnames(this_mat) %in% all_comparisons])
        tmp <- tmp[,match(all_comparisons, colnames(tmp))]
	colnames(tmp) <- all_comparisons
	tmp <- cbind(this_mat[,1:11], tmp)

        present <- !( colSums(is.na(tmp)) == nrow(tmp))

	mat_sum <- convert_cellphonedbe_2_local(tmp, expr_threshold=-1)
	mat_n <- convert_cellphonedbe_2_local(tmp, expr_threshold=0.5)
	if (is.null(dim(sum_exprLevel)) ) {
		sum_exprLevel <- mat_sum
		sum_interPairs <- mat_n
		sum_present <- present
	} else {
		sum_exprLevel <- sum_exprLevel+mat_sum
		sum_interPairs <- sum_interPairs+mat_n
        	sum_present <- sum_present + present
	}
}
types <- get_cell_types(names(sum_present))

mat_present <- 0*sum_exprLevel;
for(r in 1:nrow(types)) {
	mat_present[rownames(sum_exprLevel) == types[r,1], colnames(sum_exprLevel) == types[r,2]] <- sum_present[r]
}
avg_exprLevel <- sum_exprLevel/mat_present
avg_interPairs <- sum_interPairs/mat_present

# low-limit
avg_interPairs[avg_interPairs<3] <- 0;
avg_exprLevel[avg_interPairs==0] <- 0;

# cell_type_specific

require(igraph)
graph <- graph_from_adjacency_matrix(avg_exprLevel, mode="undirected", weighted=TRUE, diag=TRUE)
edge_colour <- as.vector(avg_interPairs);
edge_colour <- edge_colour[edge_colour > 0];
greys <- colorRampPalette(c("grey85", "black"))(10)
bins <- seq(from=range(edge_colour)[1], to=range(edge_colour)[2], length=11)
bins[1] <- bins[1]-1
bins[11] <- bins[11]+1
edge_colour_vec <- greys[cut(edge_colour, bins)]

source("~/scripts/LiverMap2.0//Colour_Scheme.R")
rownames(Cell_type_colours) <- Cell_type_colours[,1]

vertex_colour_vec <- c(rep(Cell_type_colours["AntiBcell",2], 6), 
			Cell_type_colours["Hepatocyte",2], 
			rep(Cell_type_colours["Cholangiocyte",2], 3), 
			rep(Cell_type_colours["LSECs",2], 5),
			rep(Cell_type_colours["Hepatocyte",2], 9), 
			rep(Cell_type_colours["InfMac",2],2), 
			rep(Cell_type_colours["NonInfMac",2], 6), 
			rep(Cell_type_colours["Tcell",2], 2),
			Cell_type_colours["NKcells",2],
			rep(Cell_type_colours["Tcell",2], 3),
			Cell_type_colours["NKcells",2],
			Cell_type_colours["AntiBcell",2],
			rep(Cell_type_colours["Tcell",2], 2),
			rep(Cell_type_colours["Hepatocyte",2], 2),
			rep(Cell_type_colours["Stellate",2], 5)
			)

set.seed(2901)
png("Cell_cell_graph.png", width=8, height=8, units="in", res=300)
plot(graph, edge.color=edge_colour_vec, edge.width=E(graph)$weight, vertex.color=vertex_colour_vec, layout=layout.circle)
dev.off()

saveRDS(list(avg_exprLevel = avg_exprLevel, avg_interPairs=avg_interPairs, vertex.color=vertex_colour_vec), file="Cell_cell_graph.rds")

# Graph 2 = Ligand -> Receptor
sum_exprLevel <- c();
sum_interPairs  <- c();
sum_present <- c(); # get this from n_presnt

for (i in 1:length(all_mats)) {
        this_mat <- all_mats[[i]]
        this_mat <- this_mat[match( all_interactors[,1], this_mat[,1]),]
        rownames(this_mat) <- all_interactors[,1]
        tmp <- as.matrix(this_mat[,colnames(this_mat) %in% all_comparisons])
        tmp <- tmp[,match(all_comparisons, colnames(tmp))]
        colnames(tmp) <- all_comparisons
        tmp <- cbind(this_mat[,1:11], tmp)

        present <- !( colSums(is.na(tmp)) == nrow(tmp))

        mat_sum <- convert_cellphonedbe_2_ligand_receptor(tmp, expr_threshold=-1)
        mat_n <- convert_cellphonedbe_2_ligand_receptor(tmp, expr_threshold=0.5)
        if (is.null(dim(sum_exprLevel)) ) {
                sum_exprLevel <- mat_sum
                sum_interPairs <- mat_n
                sum_present <- present
        } else {
                sum_exprLevel <- sum_exprLevel+mat_sum
                sum_interPairs <- sum_interPairs+mat_n
                sum_present <- sum_present + present
        }
}
types <- get_cell_types(names(sum_present))

mat_present <- 0*sum_exprLevel;
for(r in 1:nrow(types)) {
        mat_present[rownames(sum_exprLevel) == types[r,1], colnames(sum_exprLevel) == types[r,2]] <- sum_present[r]
}
avg_exprLevel <- sum_exprLevel/mat_present
avg_interPairs <- sum_interPairs/mat_present

# low-limit
avg_interPairs[avg_interPairs<1.5] <- 0;
avg_exprLevel[avg_interPairs==0] <- 0;

require(igraph)
graph <- graph_from_adjacency_matrix(avg_exprLevel, mode="directed", weighted=TRUE, diag=TRUE)
edge_colour <- as.vector(avg_interPairs);
edge_colour <- edge_colour[edge_colour > 0];
greys <- colorRampPalette(c("grey85", "black"))(10)
bins <- seq(from=range(edge_colour)[1], to=range(edge_colour)[2], length=11)
bins[1] <- bins[1]-1
bins[11] <- bins[11]+1
edge_colour_vec <- greys[cut(edge_colour, bins)]

source("~/scripts/LiverMap2.0//Colour_Scheme.R")
rownames(Cell_type_colours) <- Cell_type_colours[,1]

vertex_colour_vec <- c(rep(Cell_type_colours["AntiBcell",2], 6),
                        Cell_type_colours["Hepatocyte",2],
                        rep(Cell_type_colours["Cholangiocyte",2], 3),
                        rep(Cell_type_colours["LSECs",2], 5),
                        rep(Cell_type_colours["Hepatocyte",2], 9),
                        rep(Cell_type_colours["InfMac",2],2),
                        rep(Cell_type_colours["NonInfMac",2], 6),
                        rep(Cell_type_colours["Tcell",2], 2),
                        Cell_type_colours["NKcells",2],
                        rep(Cell_type_colours["Tcell",2], 3),
                        Cell_type_colours["NKcells",2],
                        Cell_type_colours["AntiBcell",2],
                        rep(Cell_type_colours["Tcell",2], 2),
                        rep(Cell_type_colours["Hepatocyte",2], 2),
                        rep(Cell_type_colours["Stellate",2], 5)
                        )

set.seed(2901)
png("Cell_cell_graph_lig_receptor.png", width=8, height=8, units="in", res=300)
plot(graph, edge.color=edge_colour_vec, edge.width=E(graph)$weight, vertex.color=vertex_colour_vec, layout=layout.circle)
dev.off()
saveRDS(list(avg_exprLevel = avg_exprLevel, avg_interPairs=avg_interPairs, vertex.color=vertex_colour_vec), file="Cell_cell_lig_receptor.rds")


## 1 cell type vs others plots
expr_threshold = 1; # expression level required of the interaction to keep it

# Boxplot of observed - Expected interacting pairs across all samples by partner cell-types.
