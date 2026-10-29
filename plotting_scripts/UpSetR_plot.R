#Author: Lars Moret, MSc student @ FunGen
#Date: 29/10/2026
#Description:
#    R script for visualizing clonal-groups and ther CoO (country of origin)
#Inputs: afumigatus_clonal.genome and the metadata file
## !! Please note that this is an adapted script from Snelders and colleagues, and I did not write this by myself !!

## load libraries
library(dplyr)
library(igraph)
library(readxl)
library(UpSetR)

## load files
global_diversity <- read.table("C:/Users/larsm/Downloads/afumigatus_clonal.genome", header = TRUE)
global_metadata <- read_excel("C:/Users/larsm/Downloads/Papers_included_metadata.xlsx") %>%
  mutate(sample_code = trimws(as.character(sample_code)), country_origin = trimws(as.character(country_origin)), source = trimws(as.character(source)))

clonal_data <- global_diversity %>% filter(DST >= 0.998) ## This sets the distance cutoff from PLINK, in this case the cutoff is both based on this study and Snelders et al, 2025
g <- graph_from_data_frame(clonal_data %>% select(IID1, IID2), directed = FALSE)

components_result <- components(g)
clonal_isolates <- data.frame(sample_id = names(components_result$membership), component_id = as.character(components_result$membership), stringsAsFactors = FALSE)
group_sizes <- table(clonal_isolates$component_id)
ordered_components <- names(sort(group_sizes, decreasing = TRUE))

group_mapping <- setNames(sprintf("ClonalGroup_%03d", seq_along(ordered_components)), ordered_components)

clonal_isolates <- clonal_isolates %>%
  mutate(clonal_group = unname(group_mapping[component_id]), group_size = as.integer(group_sizes[component_id])) %>% select(sample_id, clonal_group, group_size)

cat(n_distinct(clonal_isolates$clonal_group), "total clonal groups")
cat("With", nrow(clonal_isolates), "Isolates in these clonal groups")

## summarize distributions
group_size_summary <- clonal_isolates %>% distinct(clonal_group, group_size) %>% arrange(desc(group_size))
#write.csv(group_size_summary, "clonal_group_size_summary.csv", row.names = FALSE)

# Extract the SRA/DRA/ENA run accession from each isolate name and match names to metadata
clonal_metadata <- clonal_isolates %>%
  mutate(accession_id = ifelse(grepl("(SRR|ERR|DRR)[0-9]+", sample_id),
      sub(".*?((?:SRR|ERR|DRR)[0-9]+).*", "\\1", sample_id, perl = TRUE), NA_character_)) %>%
  left_join(global_metadata %>%
      transmute(accession_id = sample_code, country = country_origin, source = source) %>%
      filter(!is.na(accession_id), accession_id != "") %>%
      distinct(accession_id, .keep_all = TRUE), by = "accession_id") %>%
  mutate(country_lower = tolower(trimws(country)), source_lower = tolower(trimws(source)))

cat(sum(!is.na(clonal_metadata$country)), "isolates matched to metadata")
cat(sum(is.na(clonal_metadata$country) | clonal_metadata$country == ""), "isolats without metadata")


clonal_metadata <- clonal_metadata %>%
  mutate(country_set = case_when(
      country_lower == "netherlands" ~ "Netherlands",
      country_lower %in% c(
        "uk", "united kingdom", "great britain",
        "england", "scotland", "wales",
        "northern ireland") ~ "UK",
      
      country_lower %in% c(
        "usa", "us", "united states",
        "united states of america") | grepl("^usa\\s*\\(", country_lower) ~ "USA",
      
      country_lower == "germany" ~ "Germany",
      country_lower == "japan" ~ "Japan",
      country_lower == "republic of ireland" ~ "Ireland",
      country_lower == "new zealand" ~ "New Zealand",
      country_lower == "iceland" ~ "Iceland",
      country_lower == "india" ~ "India",
      country_lower == "australia" ~ "Australia",
      country_lower == "belgium" ~ "Belgium",
      country_lower == "china" ~ "China",
      country_lower == "colombia" ~ "Colombia",
      country_lower == "thailand" ~ "Thailand",
      country_lower == "zimbabwe" ~ "Zimbabwe",
      country_lower == "cuba" ~ "Cuba",
      country_lower == "mauritius" ~ "Mauritius",
      country_lower == "portugal" ~ "Portugal",
      country_lower == "russia" ~ "Russia",
      country_lower == "taiwan" ~ "Taiwan",
      TRUE ~ NA_character_),
    source_set = case_when(
      grepl("^clinical($|;)", source_lower) ~ "clinical",
      source_lower == "environmental" ~ "environmental",
      source_lower == "animal" ~ "animal",
      TRUE ~ NA_character_))

## define country and souce sets, a group belongs to a set if 1 or more isolate in it matches that country/source combo
set_definitions <- list(
  Netherlands_clinical = c("Netherlands", "clinical"),
  Netherlands_environmental = c("Netherlands", "environmental"),
  UK_clinical = c("UK", "clinical"),
  UK_environmental = c("UK", "environmental"),
  USA_clinical = c("USA", "clinical"),
  USA_environmental = c("USA", "environmental"),
  Germany_clinical = c("Germany", "clinical"),
  Germany_environmental = c("Germany", "environmental"),
  Japan_clinical = c("Japan", "clinical"),
  Ireland_clinical = c("Ireland", "clinical"),
  Ireland_environmental = c("Ireland", "environmental"),
  NewZealand_clinical = c("New Zealand", "clinical"),
  NewZealand_environmental = c("New Zealand", "environmental"))

listInput <- lapply(set_definitions, function(category) { clonal_metadata %>%
    filter(country_set == category[1], source_set == category[2], !is.na(clonal_group)) %>%
    pull(clonal_group) %>%
    unique()})

## compute set sizes and remove empty catagories
set_sizes <- data.frame(Set = names(listInput), Number_of_clonal_groups = lengths(listInput))
#print(set_sizes)
write.csv(set_sizes ,"clonal_group_upset_set_sizes.csv", row.names = FALSE)
listInput <- listInput[lengths(listInput) > 0]
display_names <- gsub("_", " ", names(listInput))
display_names <- gsub("NewZealand", "New Zealand", display_names)
display_names <- paste0(display_names," (n = ", lengths(listInput),")")
names(listInput) <- display_names
print(display_names)

## create the upsetR plot from Snelder et al 2025
upset_data <- UpSetR::fromList(listInput)
pdf("clonal_group_upset_plot.pdf", width = 14, height = 7)
UpSetR::upset(
  upset_data,
  order.by = "freq",
  nsets = length(listInput),
  nintersects = NA,
  keep.order = TRUE,
  sets = names(listInput),
  mainbar.y.label = "Number of clonal groups",
  sets.x.label = "",
  sets.bar.color = "white",
  show.numbers = "no",
  text.scale = 1.1)

dev.off()
