#Author: Lars Moret
#Description: Plot the "flavours" of MSH6, as well as the total SNP's per flavour
#Usage: - input: an annotated GFF file, 
#       - the gene variant matrix (produced by gene_variant_matrix_v2.py)
#       - the consequence file (produced by final_gene_consequences.py)
#       - the total snps file (produced by total_snps.py)

## Load libraries
#install.packages("ggplot2")

## For the gene_name, gene_id and consequence file, look at this lower part:
#gene	gene_name	
#PMS1	AFUA_2G13410	
#MSH2	AFUA_3G09850	
#MSH3	AFUA_7G04480	
#MSH6	AFUA_4G08300
#mreA	AFUA_6G03620

library(ggplot2)

setwd("D:/GEN_matrices/rfiles")
gene_id <- "AFUA_2G13410"
gene_name <- "PMS1"
variant_matrix_file <- "gene_variant_matrix.tsv"
consequence_file_name <- "PMS1_consequences.tsv"
snp_file_name <- "total_snps.tsv"
flavour_threshold <- 10 ## controls the amount of samples within a flavour to pass
gff_file <- read.table("genes.gff", header = FALSE,sep = "\t",comment.char = "#",quote = "",stringsAsFactors = FALSE)
colnames(gff_file) <- c("chrom","source","type","start","end","score","strand","phase","attributes")
gene_variant_matrix <- read.table(variant_matrix_file,header = TRUE,row.names = 1,sep = "\t",check.names = FALSE)
consequence_file <- read.table(consequence_file_name,header = TRUE,sep = "\t",stringsAsFactors = FALSE,quote = "",check.names = FALSE)
snp_counts <- read.table(snp_file_name,header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)

gene <- gff_file[gff_file$type == "gene" &grepl(gene_id,gff_file$attributes,fixed = TRUE),]
gene_chrom <- gene$chrom[1]
gene_start <- gene$start[1]
gene_stop <- gene$end[1]
gene_strand <- gene$strand[1]

gene_exons <- gff_file[gff_file$type == "exon" & grepl(gene_id,gff_file$attributes,fixed = TRUE),]
gene_exons <- gene_exons[order(gene_exons$start),]
exons_plot <- data.frame(start = gene_exons$start, end = gene_exons$end)

gene_cds <- gff_file[gff_file$type == "CDS" & grepl(gene_id,gff_file$attributes,fixed = TRUE),]
gene_cds <- gene_cds[order(gene_cds$start),]


if (nrow(gene_exons) > 1) {
  introns_plot <- data.frame(start = head(gene_exons$end, -1), end = tail(gene_exons$start, -1))
} else {
  introns_plot <- data.frame(start = numeric(0), end = numeric(0))
}


if (nrow(gene_cds) > 0) {
  cds_start <- min(gene_cds$start)
  cds_end <- max(gene_cds$end)
  
  if (gene_strand == "+") {
    utr_plot <- data.frame(type = c("5' UTR","3' UTR"),start = c(gene_start,cds_end + 1),end = c(cds_start - 1,gene_stop))
  } else {
    utr_plot <- data.frame(type = c("5' UTR","3' UTR"),start = c(cds_end + 1,gene_start),end = c(gene_stop,cds_start - 1))
  }
  utr_plot <- utr_plot[utr_plot$start <= utr_plot$end,]
  
} else {
  utr_plot <- data.frame(type = character(0),start = numeric(0),end = numeric(0))
}

## plot 500 bp before and after the gene
plot_start <- gene_start - 500
plot_end <- gene_stop + 500

variant_names <- colnames(gene_variant_matrix)
variant_data <- data.frame(variant = variant_names,stringsAsFactors = FALSE)

variant_parts <- strsplit(variant_data$variant,":",fixed = TRUE)
variant_data$chrom <- sapply(variant_parts, function(x) {
    if (length(x) >= 3) {x[2]
    } else {NA}})

variant_data$position_string <- sapply(variant_parts, function(x) {
    if (length(x) >= 3) {x[3]
    } else {NA}})


variant_data$P <- as.numeric(sub("_.*","",variant_data$position_string))
variant_data$change <- sapply(variant_parts,function(x) {
    if (length(x) >= 3) {sub("^[0-9]+_","",x[3])
    } else {NA}})

number_on_gene_chrom <- sum(variant_data$chrom == gene_chrom,na.rm = TRUE)

variant_data$region <- "intergenic"
correct_chromosome <- variant_data$chrom == gene_chrom

variant_data$region[correct_chromosome & variant_data$P >= gene_start & variant_data$P <= gene_stop] <- "intron"

if (nrow(utr_plot) > 0) {for (i in seq_len(nrow(utr_plot))) {
    variant_data$region[correct_chromosome & variant_data$P >= utr_plot$start[i] & variant_data$P <= utr_plot$end[i]] <- utr_plot$type[i]}}

if (nrow(gene_cds) > 0) {is_cds <- sapply(variant_data$P,function(p) {
      any(p >= gene_cds$start & p <= gene_cds$end)})

  variant_data$region[correct_chromosome & is_cds] <- "CDS"}


flavours <- aggregate(rep(1, nrow(gene_variant_matrix)), by = as.data.frame(gene_variant_matrix), FUN = sum)
names(flavours)[ncol(flavours)] <- "Number_of_isolates"

flavours$Flavour <- paste0("Flavour_",seq_len(nrow(flavours)))
flavours_plot <- flavours[flavours$Number_of_isolates > flavour_threshold,]
flavours_plot <- flavours_plot[order(flavours_plot$Number_of_isolates,decreasing = TRUE),]
flavour_long <- reshape(flavours_plot,varying = variant_names,v.names = "present",timevar = "variant",times = variant_names,idvar = c("Flavour","Number_of_isolates"),direction = "long")
flavour_long <- flavour_long[flavour_long$present == 1,]
flavour_long$P <- variant_data$P[match(flavour_long$variant,variant_data$variant)]
flavour_long$chrom <- variant_data$chrom[match(flavour_long$variant,variant_data$variant)]
flavour_long$region <- variant_data$region[match(flavour_long$variant,variant_data$variant)]
flavour_long <- flavour_long[flavour_long$chrom == gene_chrom &flavour_long$P >= plot_start &flavour_long$P <= plot_end,]
flavour_levels <- flavours_plot$Flavour
flavour_long$Flavour <- factor(flavour_long$Flavour,levels = rev(flavour_levels))
flavour_long$Flavour <- droplevels(flavour_long$Flavour)
flavour_long$y <- as.numeric(flavour_long$Flavour)
flavour_long$variant_match <- sub("^[^:]+:","",flavour_long$variant)

consequence_file$variant_match <- sub("^[^:]+:","",consequence_file$variant)
track_list <- split(flavour_long$P,flavour_long$Flavour)
flavour_tracks <- data.frame(Flavour = names(track_list),start = sapply(track_list,min),end = sapply(track_list,max),row.names = NULL)
flavour_tracks$y <- match(flavour_tracks$Flavour,levels(flavour_long$Flavour))
consequence_match <- match(flavour_long$variant_match,consequence_file$variant_match)
flavour_long$consequence <- consequence_file$consequence[consequence_match]
flavour_long$hgvs_c <- consequence_file$hgvs_c[consequence_match]
flavour_long$hgvs_p <- consequence_file$hgvs_p[consequence_match]


## Specify the variant classes (e.g. synonymous, non-synonymous, etc.) from snpEff
flavour_long$variant_class <- "other"
flavour_long$variant_class[grepl("synonymous_variant",flavour_long$consequence)] <- "synonymous"
flavour_long$variant_class[grepl("missense_variant",flavour_long$consequence)] <- "missense"
flavour_long$variant_class[grepl("frameshift_variant",flavour_long$consequence)] <- "frameshift"
flavour_long$variant_class[grepl("stop_gained|stop_lost|start_lost",flavour_long$consequence)] <- "nonsense"
flavour_long$variant_class[grepl("inframe_insertion|inframe_deletion",flavour_long$consequence)] <- "inframe_indel"
flavour_long$variant_class[grepl("intron_variant",flavour_long$consequence)] <- "intronic"
flavour_long$variant_class[grepl("5_prime_UTR_variant",flavour_long$consequence)] <- "5_UTR"
flavour_long$variant_class[grepl("3_prime_UTR_variant",flavour_long$consequence)] <- "3_UTR"

all_variants <- unique(flavour_long[,c("variant","P","region","consequence","hgvs_c","hgvs_p")])
all_variants <- all_variants[order(all_variants$P),]

flavours$Flavour <- paste0("Flavour_", seq_len(nrow(flavours)))

## manual labor; rename these to change the "flavours" in the plot to a more presentable name
flavour_names <- c(
  "Flavour_421"  = "WT",
  "Flavour_319" = "PMS1_variant_1",
  "Flavour_100" = "PMS1_variant_2",
  "Flavour_242" = "PMS1_variant_3")

flavours$Flavour <- ifelse(flavours$Flavour %in% names(flavour_names),flavour_names[flavours$Flavour],flavours$Flavour)

## Plot the gene with its variations
gene_plot <- ggplot() + geom_segment(aes(x = plot_start,xend = plot_end,y = 0,yend = 0),linewidth = 0.6,color = "grey50") +
  geom_rect(data = introns_plot,aes(xmin = start,xmax = end,ymin = -0.15,ymax = 0.15),fill = "steelblue3",color = NA) +
  geom_rect(data = exons_plot,aes(xmin = start,xmax = end,ymin = -0.25,ymax = 0.25),fill = "steelblue3",color = NA) +
  geom_rect(data = utr_plot,aes(xmin = start,xmax = end,ymin = -0.25,ymax = 0.25),fill = "orange",color = NA) +
  geom_segment(data = flavour_tracks,aes(x = plot_start,xend = plot_end,y = -y,yend = -y),linewidth = 0.5,color = "black") +
  geom_segment(data = flavour_long,aes(x = P,xend = P,y = -y + 0.12,yend = -y - 0.12,color = variant_class),linewidth = 0.8) +
  scale_color_manual(values = c(synonymous = "green3", missense = "red",frameshift = "purple",nonsense = "black",inframe_indel = "blue",intronic = "orange",`5_UTR` = "gold",`3_UTR` = "goldenrod",other = "grey50"),
  breaks = c("synonymous","missense","frameshift","nonsense","inframe_indel","intronic","5_UTR","3_UTR","other"),
  labels = c("Synonymous","Non-synonymous","Frameshift","Stop/start","In-frame indel","Intronic","5′ UTR","3′ UTR","Other")) +
  geom_text(data = flavour_tracks,aes(x = plot_start - 25,y = -y,label = Flavour),hjust = 1,size = 3) +
  coord_cartesian(xlim = c(plot_start,plot_end),ylim = c(-max(flavour_long$y) - 0.5,0.5)) +
  labs(title = paste(gene_name,"flavours"),x = paste("Genomic position (",gene_chrom,")"),y = NULL,color = "Variant consequence") +
  theme_classic() + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())

print(gene_plot)

flavours <- aggregate(rep(1,nrow(gene_variant_matrix)),by = as.data.frame(gene_variant_matrix),FUN = sum)
names(flavours)[ncol(flavours)] <- "Number_of_isolates"
flavours$Flavour <- paste0("Flavour_",seq_len(nrow(flavours)))
isolate_patterns <- apply(gene_variant_matrix,1,paste,collapse = "_")
flavour_patterns <- apply(flavours[,variant_names,drop = FALSE],1,paste,collapse = "_")
isolate_flavours <- data.frame(isolate = rownames(gene_variant_matrix),pattern = isolate_patterns,stringsAsFactors = FALSE)
isolate_flavours$Flavour <- flavours$Flavour[match(isolate_flavours$pattern,flavour_patterns)]
flavour_sizes <- table(isolate_flavours$Flavour)
isolate_flavours_plot <- isolate_flavours[isolate_flavours$Flavour %in% names(flavour_sizes[flavour_sizes > flavour_threshold]),]

snp_flavour <- merge(snp_counts,isolate_flavours_plot[,c("isolate","Flavour")],by = "isolate")

## plot the per "flavour" total snps (note this is not filtered to clonal groups, these are raw population plots and not account for diversity)
snp_boxplot <- ggplot(snp_flavour,aes(x = Flavour,y = total_SNPs)) +
  geom_boxplot() +labs(title = paste("Total SNPs per isolate by", gene_name, "flavour"),
    x = paste(gene_name,"flavour"),
    y = "Total SNPs per isolate") + theme_classic()

print(snp_boxplot)
