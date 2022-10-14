library(stringr)
library(dplyr)
library(ggplot2)
library(igraph)
library(tidyr)

# Number of ISs
info_files <- list.files("data/ISC/protocols/palidis_output/v3.1.0", pattern = "*_insertion_sequences_info.txt", full.names = TRUE)
is_info <- data.frame()
for (i in 1:length(info_files)) {
  info = file.info(info_files[i])
  if (info$size != 0) {
    tmp <- read.delim(info_files[i], header = TRUE, stringsAsFactors = FALSE)
    is_info <- rbind(is_info, tmp)
  }
}

length(unique(is_info$contig)) # Number of contigs with an IS
length(unique(is_info$sample_id)) # Number of samples with an IS
length(is_info$IS_name) # Number of ISs found 

write.csv(is_info, "supplementary/Supplementary_Data_2.csv", quote = FALSE, row.names = FALSE)

# Get unique ISs
isc <- read.delim("data/ISC/insertion_sequence_catalogue_info.txt", stringsAsFactors = FALSE)
sum(unique(is_info$IS_name) %in% isc$IS_name)
isc_current <- isc[isc$IS_name %in% unique(is_info$IS_name),]
isc_current$length <- as.integer(sapply(gsub("IS_length_", "", isc_current$IS_name), function(x) str_split(x, "-")[[1]][1]))

# Plot length of ISs
tiff("figures/is_length.tiff", width=2000, height=2000, res = 400)
ggplot(data = isc_current, aes(1:nrow(isc_current), sort(length))) +
  geom_point() +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, max(isc_current$length)+500, by = 500)) +
  scale_x_continuous(breaks = seq(0, nrow(isc_current)+100, by = 100)) +
  ylab("IS length") + xlab("n of ISC")
dev.off()
min(isc_current$length) # Minimum IS length
max(isc_current$length) # Maximum IS length

# Plot transposases
isc_transposase <- isc_current %>%
  mutate(description = strsplit(as.character(description), ";")) %>% 
  unnest(description) %>%
  group_by(description) %>%
  summarise(n = n())
isc_transposase$description = factor(isc_transposase$description, levels=unique(sort(isc_transposase$description, decreasing = TRUE)))
length(unique(isc_transposase$description)) # Unique transposases

tiff("figures/transposase.tiff", width=4000, height=4000, res = 400)
ggplot(data = isc_transposase, aes(description, n)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab("Number") + xlab("Transposase") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()
dev.off()

# How many are in ISfinder?
is_finder_out <- read.delim("data/isfinder_081022_blastn_out_tabular_evalue_0.01.txt", header = FALSE, stringsAsFactors = FALSE)
names(is_finder_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

categorise_blast <- function(x) {
  category <- 'None'
  if (min(x$evalue) < 0.01) {
    category <- 'Loose'
  }
  if (min(x$evalue) < 1e-50) {
    category <- 'Strict'
  }
  return(category)
}

is_finder_hits <- is_finder_out %>%
  group_by(qseqid, sseqid) %>%
  summarise(hit = categorise_blast(.data))

length(unique(is_finder_hits$qseqid[is_finder_hits$hit == "Strict"])) # Number ISs found in ISfinder (strict homology)
length(unique(is_finder_hits$qseqid[is_finder_hits$hit == "Loose"])) # Number ISs found in ISfinder (loose homology)

# ISs found in 661k pathogen genomes
ena_metadata <- read.delim("data/File9_all_metadata_ena_661K.txt", header = TRUE, stringsAsFactors = FALSE)
cobs_out <- read.delim("data/insertion_sequence_catalogue.fasta_0.9_results_table.txt", header = TRUE, stringsAsFactors = FALSE)

cobs_df <- left_join(cobs_out, ena_metadata, by = "sample_id")

length(unique(cobs_df$query)) # Number of ISs found in ENA cobs index
cobs_df <- cobs_df[!is.na(cobs_df$scientific_name),]
length(unique(cobs_df$sample_id)) # Number of unique sources from ENA cobs index

# Cleaning
cobs_df <- cobs_df[cobs_df$scientific_name != "unidentified bacterium",]

cobs_df$scientific_name <- gsub('\\[', '', cobs_df$scientific_name)
cobs_df$scientific_name <- gsub('\\] ', '', cobs_df$scientific_name)

cobs_df$scientific_name <- gsub("Clostridiuminnocuum", "Clostridium innocuum", cobs_df$scientific_name)
cobs_df$scientific_name <- gsub("Clostridiumspiroforme", "Clostridium spiroforme", cobs_df$scientific_name)
cobs_df$scientific_name <- gsub("Clostridiumsymbiosum", "Clostridium symbiosum", cobs_df$scientific_name)
cobs_df$scientific_name <- gsub("Eubacteriumrectale", "Eubacterium rectale", cobs_df$scientific_name)

cobs_df$cultured <- 'Y'
cobs_df$cultured[grep("uncultured", cobs_df$scientific_name)] <- 'N'
cobs_df$scientific_name <- gsub('uncultured ', '', cobs_df$scientific_name)
cobs_df$genus <- sapply(cobs_df$scientific_name, function(x) str_split(x, " ")[[1]][1])
cobs_df$species <- sapply(cobs_df$scientific_name, function(x) paste(str_split(x, " ")[[1]][1], str_split(x, " ")[[1]][2]))
cobs_df$species[is.na(cobs_df$species) & !is.na(cobs_df$scientific_name) & cobs_df$scientific_name != "bacterium"] <- 'sp.'

# Number of sources being isolate vs community
length(unique(cobs_df$sample_id[cobs_df$genus == "bacterium"])) # community
length(unique(cobs_df$sample_id[cobs_df$genus != "bacterium"])) # isolate

length(unique(cobs_df$genus[cobs_df$genus != "bacterium" & !grepl(" bacterium", cobs_df$species)])) # Number of unique genera
length(unique(cobs_df$species[!grepl("sp.", cobs_df$species) & !grepl("NA", cobs_df$species) & !grepl(" bacterium", cobs_df$species)])) # Number of unique known species

# Summary of genera
genera_summary <- cobs_df %>%
  filter(genus != "bacterium" & !grepl(" bacterium", species)) %>%
  group_by(genus) %>%
  summarise(n = n_distinct(sample_id))

tiff("figures/genera_sample_summary.tiff", width=1000, height=500, res = 100)
ggplot(genera_summary, aes(genus, n)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  ylab("No. samples") + xlab("Genus") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Summary of species
species_summary <- cobs_df %>%
  group_by(species) %>%
  summarise(n = n_distinct(sample_id)) %>%
  filter(species != "bacterium NA" & !grepl(" bacterium", species)) %>%
  mutate(species = gsub("NA", "sp.", species))

tiff("figures/species_sample_summary.tiff", width=2000, height=500, res = 100)
ggplot(species_summary[species_summary$species != "bacterium NA",], aes(species, n)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  ylab("No. samples") + xlab("Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Shared ISs
is_shared_genera <- cobs_df %>%
  filter(genus != "bacterium" & !grepl(" bacterium", species)) %>%
  group_by(query) %>%
  summarise(n = n_distinct(genus)) %>%
  filter(n > 1) %>%
  arrange(desc(n)) %>%
  mutate(in_isfinder = query %in% is_finder_out$qseqid)

length(is_shared_genera$query) # Number of ISs shared across genera

no_is_shared_genera <- is_shared_genera %>%
  group_by(n, in_isfinder) %>%
  summarise(n_is = n())

tiff("figures/shared_is.tiff", width=500, height=500, res = 100)
ggplot(no_is_shared_genera, aes(n, n_is, fill=in_isfinder)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  ylab("No. ISs shared across genera") + xlab("No. genera") +
  scale_x_continuous(breaks = seq(0, max(no_is_shared_genera$n), by = 1))
dev.off()

is_shared_species <- cobs_df %>%
  filter(genus != "bacterium" & !grepl(" bacterium", species) & !grepl("sp.", genus) & !grepl("NA", genus)) %>%
  group_by(query) %>%
  summarise(n = n_distinct(species)) %>%
  filter(n > 1) %>%
  arrange(desc(n))

length(is_shared_species$query) # Number of ISs shared across species

shared_most <- is_shared_genera$query[is_shared_genera$n == max(is_shared_genera$n)]
shared_most # IS that is shared most
max(is_shared_genera$n) # Number of genera
max(is_shared_species$n) # Number of species
shared_most %in% is_finder_out$qseqid # In ISfinder?
is_finder_most_shared <- is_finder_out[is_finder_out$qseqid == shared_most,] %>%
  filter(evalue == min(evalue)) %>%
  filter(pident == max(pident))
is_finder_most_shared$sseqid # IS from ISfinder

shared_most_not_isfinder <- is_shared_genera$query[!is_shared_genera$in_isfinder][is_shared_genera$n[!is_shared_genera$in_isfinder] == max(is_shared_genera$n[!is_shared_genera$in_isfinder])]
shared_most_not_isfinder # IS shared most not in ISfinder

# Genera/IS network
is_shared_named_genera <- left_join(is_shared_genera, cobs_df, by = "query") %>%
  filter(genus != "bacterium" & !grepl(" bacterium", species)) %>%
  select(query, genus)

g <- graph.data.frame(d = is_shared_named_genera, directed = FALSE)
lay <- layout_with_fr(g)
vertex_names <- V(g)$name

catgry = c("IS also in ISfinder", "IS not in ISfinder")
f <- factor(rep(catgry, length = length(V(g))))
colrs = c("#619CFF", "#F8766D")
vcols <- colrs[f]
vertex_colours <- rep("lightblue", length(vertex_names))
vertex_colours[grepl("IS_", vertex_names) & vertex_names %in% is_finder_out$qseqid] <- colrs[catgry == "IS also in ISfinder"]
vertex_colours[grepl("IS_", vertex_names) & !(vertex_names %in% is_finder_out$qseqid)] <- colrs[catgry == "IS not in ISfinder"]

vertex_shapes <- rep("circle", length(vertex_names))
vertex_shapes[vertex_names %in% c(shared_most_not_isfinder, shared_most)] <- "square"

vertex_names[grep("IS_", vertex_names)] <- ""

tiff("figures/genera_is_network.tiff", width=2000, height=2000, res = 200)
set.seed(1492)
l <- layout.fruchterman.reingold(g, niter=5000)
plot(g, layout=l, 
     vertex.size=1.5, 
     vertex.label=vertex_names, 
     vertex.color = vertex_colours, 
     vertex.label.color = "black",
     vertex.label.cex = 0.55,
     vertex.label.family = "Helvetica",
     vertex.frame.color =  NA,
     edge.width = 0.5,
     vertex.shape = vertex_shapes)
legend("bottomright", legend = levels(f), col = vcols, pch = 16, bty = "n")
dev.off()

# ISfinder representation
isfinder_summary <- cobs_df %>%
  filter(genus != "bacterium") %>%
  group_by(query) %>%
  select(query, species) %>%
  filter(species != "bacterium NA" & !grepl(" bacterium", species)) %>%
  mutate(species = gsub("NA", "sp.", species)) %>%
  mutate(in_isfinder = query %in% is_finder_out$qseqid) %>%
  group_by(species, in_isfinder) %>%
  summarise(n = n())
isfinder_summary$species = factor(isfinder_summary$species, levels=unique(sort(isfinder_summary$species, decreasing = TRUE)))

tiff("figures/isfinder_summary.tiff", width=1000, height=1600, res = 100)
ggplot(isfinder_summary, aes(n, species, fill=in_isfinder)) +
  geom_bar(position = "dodge", stat = "identity") + 
  theme_minimal() +
  scale_fill_manual("IS in ISfinder", values = c("FALSE" = colrs[catgry == "IS not in ISfinder"], "TRUE" = colrs[catgry == "IS also in ISfinder"])) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_continuous(breaks = seq(0, max(isfinder_summary$n)+10, by = 10)) +
  xlab("No. ISs")
dev.off()
  

