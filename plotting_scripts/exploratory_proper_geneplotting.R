## this file has been depricated, but can be used as exploratory before proper_geneplotting.R
##
##Author: Lars Moret
#Description: Plot the "flavours" of MSH6, as well as the total SNP's per flavour
#Usage: - input: an annotated GFF file, 
#       - the gene variant matrix (produced by gene_variant_matrix_v2.py)
#       - the consequence file (produced by gene_consequences_v2.py)
#       - the total snps file (produced by total_snps.py)

## Load libraries
#install.packages("ggplot2")
library(ggplot2)


## Set work directory and files
setwd("D:/GEN_matrices/rfiles")
gff_file <- read.table("genes.gff", header = FALSE, sep = "\t", comment.char = "#", quote = "", stringsAsFactors = FALSE)
colnames(gff_file) <- c("chrom", "source", "type", "start",  "end", "score", "strand", "phase", "attributes")
gene_variant_matrix <- read.table("gene_variant_matrix.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
msh6 <- gff_file[gff_file$type == "gene" & grepl("AFUA_2G13410", gff_file$attributes),]
consequence_file <- read.table("PMS1_consequences.tsv",header = TRUE,sep = "\t",stringsAsFactors = FALSE,quote = "",check.names = FALSE)
snp_counts <- read.table("total_snps.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)


## Define MSH6
msh6_exons <- gff_file[gff_file$type == "exon" & grepl("AFUA_2G13410", gff_file$attributes),]
msh6_exons <- msh6_exons[order(msh6_exons$start),]
msh6_exons[, c("chrom", "start", "end", "strand")]
msh6_gene_start <- msh6$start
msh6_gene_stop <- msh6$end
msh6_gene_strand <- msh6$strand
msh6_cds <- gff_file[gff_file$type == "CDS" & grepl("AFUA_2G13410", gff_file$attributes),]
msh6_cds <- msh6_cds[order(msh6_cds$start),]
msh6_cds[, c("chrom", "start", "end", "strand")]
msh6_cds_length <- sum(msh6_cds$end - msh6_cds$start + 1)
msh6_introns <- data.frame(start = head(msh6_exons$end, -1) + 1,end = tail(msh6_exons$start, -1) - 1)
msh6_utr <- data.frame(type = c("3' UTR", "5' UTR"), start = c(msh6_gene_start,max(msh6_cds$end) + 1), end = c(min(msh6_cds$start) - 1,msh6_gene_stop))

plot_start <- msh6_gene_start - 500
plot_end <- msh6_gene_stop + 500

exons_plot <- data.frame(start = msh6_exons$start,end = msh6_exons$end)
introns_plot <- data.frame(start = head(exons_plot$end, -1),end = tail(exons_plot$start, -1))
utr_plot <- msh6_utr

variant_names <- colnames(gene_variant_matrix)
variant_data <- data.frame(variant = variant_names)
variant_data$P <- as.numeric(sub(".*:([0-9]+)_.*", "\\1", variant_data$variant))
head(variant_data, 20)

variant_data$region <- "intergenic"
variant_data$region[variant_data$P >= msh6_gene_start &variant_data$P <= msh6_gene_stop] <- "intron"
variant_data$region[variant_data$P >= msh6_utr$start[msh6_utr$type == "3' UTR"] &variant_data$P <= msh6_utr$end[msh6_utr$type == "3' UTR"]] <- "3' UTR"
variant_data$region[variant_data$P >= msh6_utr$start[msh6_utr$type == "5' UTR"] & variant_data$P <= msh6_utr$end[msh6_utr$type == "5' UTR"]] <- "5' UTR"
variant_data$region[sapply(variant_data$P, function(p) {any(p >= msh6_cds$start & p <= msh6_cds$end)})] <- "CDS"
table(variant_data$region)

flavours <- aggregate(rep(1, nrow(gene_variant_matrix)), by = as.data.frame(gene_variant_matrix), FUN=sum)
names(flavours)[ncol(flavours)] <- "Number_of_isolates"
nrow(flavours)
sum(flavours$Number_of_isolates > 5)
sort(flavours$Number_of_isolates, decreasing=T)
msh6_exons[, c("chrom", "start", "end", "strand")]

flavours$Flavour <- paste0("Flavour_", seq_len(nrow(flavours)))
flavours_plot <- flavours[flavours$Number_of_isolates > 10,] ### Change to increase the number of isolates per group threshold
flavours_plot <- flavours_plot[order(flavours_plot$Number_of_isolates, decreasing = TRUE),]
flavours_plot[, c("Flavour", "Number_of_isolates")]
flavour_long <- reshape(flavours_plot, varying = variant_names, v.names = "present", timevar = "variant", 
                        times = variant_names,idvar = c("Flavour", "Number_of_isolates"), direction = "long")
flavour_long <- flavour_long[flavour_long$present == 1,]
flavour_long$P <- variant_data$P[match(flavour_long$variant,variant_data$variant)]
head(flavour_long)

flavour_levels <- flavours_plot$Flavour
flavour_long$Flavour <- factor(flavour_long$Flavour, levels = rev(flavour_levels))
flavour_long$Flavour <- droplevels(flavour_long$Flavour)
flavour_long$y <- as.numeric(flavour_long$Flavour)
track_list <- split(flavour_long$P, flavour_long$Flavour)
flavour_tracks <- data.frame(Flavour = names(track_list), start = sapply(track_list, min), end = sapply(track_list, max), row.names = NULL)
flavour_tracks$y <- match(flavour_tracks$Flavour, levels(flavour_long$Flavour))

## Old plot (without colored tracks)
#flavour_tracks$Flavour <- rownames(flavour_tracks)
#rownames(flavour_tracks) <- NULL
#flavour_tracks$y <- as.numeric(factor(flavour_tracks$Flavour,levels = levels(flavour_long$Flavour)))

#ggplot() + geom_segment(aes(x = plot_start,xend = plot_end,y = 0,yend = 0), linewidth = 0.6,color = "grey50") +
#  geom_rect(data = introns_plot,aes(xmin = start,xmax = end,ymin = -0.15,ymax = 0.15), fill = "steelblue3",color = NA) +
#  geom_rect(data = exons_plot,aes(xmin = start,xmax = end,ymin = -0.25,ymax = 0.25),fill = "steelblue3",color = NA) +
#  geom_rect(data = utr_plot,aes(xmin = start,xmax = end,ymin = -0.25,ymax = 0.25),fill = "orange",color = NA) +
#  geom_segment(data = flavour_tracks,aes(x = start,xend = end,y = -y,yend = -y),linewidth = 0.5,color = "black") +
#  geom_segment(data = flavour_long,aes(x = P,xend = P,y = -y + 0.12,yend = -y - 0.12),linewidth = 0.8,color = "black") +
#  coord_cartesian(xlim = c(plot_start, plot_end),ylim = c(-max(flavour_long$y) - 0.5,0.5)) + labs(x = "Genomic position (NC_007197.1)",y = NULL) +
#  theme_classic() + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(), axis.title.y = element_blank())

flavour_long$consequence <- consequence_file$consequence[match(flavour_long$variant,consequence_file$variant)]
flavour_long$hgvs_c <- consequence_file$hgvs_c[match(flavour_long$variant,consequence_file$variant)]
flavour_long$hgvs_p <- consequence_file$hgvs_p[match(flavour_long$variant,consequence_file$variant)]


flavour_long$variant_class <- "other"
flavour_long$variant_class[grepl("synonymous_variant", flavour_long$consequence)] <- "synonymous"
flavour_long$variant_class[grepl("missense_variant", flavour_long$consequence)] <- "missense"
flavour_long$variant_class[grepl("frameshift_variant", flavour_long$consequence)] <- "frameshift"
flavour_long$variant_class[grepl("stop_gained|stop_lost|start_lost", flavour_long$consequence)] <- "nonsense"
flavour_long$variant_class[grepl("inframe_insertion|inframe_deletion",flavour_long$consequence)] <- "inframe_indel"
flavour_long$variant_class[grepl("intron_variant", flavour_long$consequence)] <- "intronic"
flavour_long$variant_class[grepl("5_prime_UTR_variant", flavour_long$consequence)] <- "5_UTR"
flavour_long$variant_class[grepl("3_prime_UTR_variant", flavour_long$consequence)] <- "3_UTR"


## Define the G233A variant of MSH6 (Wanted to plot this specifically)
g233a_variant <- "AFUA_4G08300:NC_007197.1:2148956_C>G"
g233a_rows <- flavour_long[flavour_long$variant == g233a_variant,,drop = FALSE]

if (nrow(g233a_rows) == 0) {
  cat("G233A was NOT found in flavour_long.\n\n")
} else {cat("G233A was found.\n\n")print(unique(g233a_rows[,c("variant","P","consequence","hgvs_c","hgvs_p")]),row.names = FALSE)}

flavour_long$plot_class <- flavour_long$variant_class
flavour_long$plot_class[flavour_long$variant == g233a_variant] <- "G233A"

## Identify all the flavours with G233A
g233a_flavours <- unique(flavour_long$Flavour[flavour_long$variant == g233a_variant])
cat("Flavours containing G233A:",length(g233a_flavours))

## And how many isolates are present in those flavours
g233a_flavour_info <- flavours_plot[flavours_plot$Flavour %in% g233a_flavours,c("Flavour", "Number_of_isolates")]
print(g233a_flavour_info,row.names = FALSE)
cat("Total isolates represented by G233A flavours:", sum(g233a_flavour_info$Number_of_isolates))

all_variants <- unique(flavour_long[,c("variant","P","consequence","hgvs_c","hgvs_p")])
all_variants$change <- sub(".*_","",all_variants$variant)
all_variants <- all_variants[order(all_variants$P),c("P","change","consequence","hgvs_c","hgvs_p")]

synonymous <- all_variants[grepl("synonymous_variant",all_variants$consequence),]
nonsynonymous <- all_variants[grepl("missense_variant|frameshift_variant|stop_gained|stop_lost|start_lost|inframe",all_variants$consequence),]

## Plot everything using GGplot2
ggplot() + geom_segment(aes(x = plot_start,xend = plot_end,y = 0,yend = 0), linewidth = 0.6,color = "grey50") +
  geom_rect(data = introns_plot,aes(xmin = start,xmax = end,ymin = -0.15,ymax = 0.15),fill = "steelblue3",color = NA) +
  geom_rect(data = exons_plot,aes(xmin = start,xmax = end,ymin = -0.25,ymax = 0.25), fill = "steelblue3",color = NA) +
  geom_rect(data = utr_plot,aes(xmin = start,xmax = end,ymin = -0.25,ymax = 0.25), fill = "orange",color = NA) +
  geom_segment(data = flavour_tracks,aes(x = start,xend = end,y = -y,yend = -y),linewidth = 0.5,color = "black") +
  geom_segment(data = flavour_long,aes(x = P,xend = P,y = -y + 0.12,yend = -y - 0.12,color = plot_class), linewidth = 0.8) +
  scale_color_manual(values = c(synonymous = "green3",missense = "red", intronic = "orange",`5_UTR` = "gold",`3_UTR` = "goldenrod", G233A = "purple"),
                     breaks = c("synonymous","missense","intronic","5_UTR","3_UTR","G233A"), 
                     labels = c("Synonymous", "Non-synonymous", "Intronic", "5′ UTR", "3′ UTR", "G233A")) +
  coord_cartesian(xlim = c(plot_start, plot_end), ylim = c(-max(flavour_long$y) - 0.5, 0.5)) +
  labs(x = "Genomic position (NC_007197.1)", y = NULL, color = "Variant consequence") + 
  geom_text(data = flavour_tracks,aes(x = plot_start - 25, y = -y, label = Flavour),hjust = 1,size = 3)
theme_classic() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())


## Total SNP per flavour comparison
## Assign the isolates in the total_snps file to their flavour
flavours <- aggregate(rep(1,nrow(gene_variant_matrix)), by= as.data.frame(gene_variant_matrix),FUN =sum)
names(flavours)[ncol(flavours)] <- "Number_of_isolates" 
flavours$Flavour <- paste0("Flavour_", seq_len(nrow(flavours)))
isolate_patterns <- apply(gene_variant_matrix, 1, paste, collapse = "_")

flavour_patterns <- apply(flavours[, variant_names, drop = FALSE],1,paste,collapse = "_")
isolate_flavours <- data.frame(isolate = rownames(gene_variant_matrix), pattern = isolate_patterns, stringsAsFactors = FALSE) 
isolate_flavours$Flavour <- flavours$Flavour[match(isolate_flavours$pattern, flavour_patterns)]
#head(isolate_flavours) 
#table(isolate_flavours$Flavour) 
#sum(is.na(isolate_flavours$Flavour))
flavour_sizes <- table(isolate_flavours$Flavour)
isolate_flavours_plot <- isolate_flavours[isolate_flavours$Flavour %in% names(flavour_sizes[flavour_sizes > 10]),]
snp_flavour <- merge(snp_counts,isolate_flavours_plot[, c("isolate", "Flavour")], by = "isolate")

ggplot(snp_flavour, aes(x = Flavour, y = total_SNPs)) +
  geom_boxplot() + labs(x = "MSH6 flavour",y = "Total SNPs per isolate") + theme_classic()

