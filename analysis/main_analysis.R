#####################
# 
#  Notes on the origin of this script
#
#####################
# This script has been adapted from a script written to handle 16S rRNA amplicon data (V3-V4 region).
# Original sample type: boar semen - but could be used for any low-microbial biomass type sample data.
#
# Originally Samples were run on MiSeq 300 pe Illumina
# In addition to semen, multiple controls were included in the original run: 
# Samples from adjacent environmental niches to act as a comparison to semen.
# I recommend the use of at least one negative control, and preferably more - e.g.: PBS, water, etc.
# Finally there are several mocks.
# 1 Mock is a logarithmic DNA standard from zymoBIOMICS.
# 2 repeats of zymoBIOMICS mock containing equal ratios of 8 organisms.
# A cell based mock created by me, in 3 dilutions.
# Undiluted contains approx 10^8 cells. Then dilutions equal to: 10^5 and 10^3 cells.
# You could also use commercial cell-based mocks.
# 
# Additional controls that could be included, but are not present in the workflow are
# spike-in controls - these are used to estimate absolute microbial quantity.
####################################################################
#
# Original pre-processing steps
#
####################################################################
# Proceeding this script, a shell-script preprocessing Qiime2 workflow was used to: 
# Trim primers and adapters with the cutadapt plugin.
# 
# Additionally, as this original run contained a high degree of off target amplification,
# Bowtie2 was used to align and remove sequences that mapped against boar genome. 
# Some samples now had ++ low sequencing depth and the script include steps to assess depth.
#
# The Qiime2-2024 amplicon workflow was then used to:
# Create a demultiplexed qza object containing read pairs.
# Truncate reads at 280f and 220r and denoise using the dada2 algorithm
# (Most of my reads were lost at quality filter, which is normal for Illumina data.)
# (Most of my reads were retained at merging and chimeric filter - all good signs.)
# Create a feature table and representative sequences again with dada2
# Create a phylogenetic tree using the mafft-tree method after first aligning all sequences
# Perform Taxonomic assignment using sklearn and pre-trained naive-bayesian classifiers.
# I used a full length SILVA 138.2 classifier, but other options are available.
#########################################################################
# this clears the global environment
rm(list = ls())
########################################################################
# Load configuration
source("config/config.R")
# Check all expected files are in place
stopifnot(
  file.exists(feature_table_qza),
  file.exists(tree_qza),
  file.exists(metadata_tsv),
  file.exists(taxonomy_qza)
)
########################################################################
# Read in the necessary packages
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(vegan)
library(patchwork)
library(grid)
library(gridExtra)
library(ade4)
library(reshape2)
library(ggrepel)
library(caret)
library(randomForest)
library(pls)
library(viridis)
library(RColorBrewer)
library(decontam)
library(ggsci)
library(microbiome)
library(ggvenn)
library(ggVennDiagram)
library(ggpubr)
library(pheatmap)
#######################################################################
# SECTION 1: Import QIIME2 data and create phyloseq object
########################################################################

# Read the qiime feature table into R
ASVs <- read_qza(feature_table_qza)
# Read the metadata file into R, ignore the filepaths to fastq files.
Metadata <- readr::read_tsv(metadata_tsv)

# Read the qiime taxonomy files into R
Taxonomy <- read_qza(taxonomy_qza)

# Parse the taxonomy string in this file to make it more manageable
Taxonomy <- parse_taxonomy(Taxonomy$data)

# Remove the d__ prefix from the Kingdom data
Taxonomy$Kingdom <- sub('^d__', '', Taxonomy$Kingdom)

# Read the tree file into R
Tree <- read_qza(tree_qza)

# Create your final Taxonomy dataframe
Taxonomy <- Taxonomy |>
  tibble::rownames_to_column("ASV_ID") |>
  select(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)
rownames(Taxonomy) <- Taxonomy$ASV_ID
Taxonomy$ASV_ID <- NULL
colnames(Taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Import into phyloseq
physeq <- qza_to_phyloseq(
  features = feature_table_qza,
  tree = tree_qza,
  metadata = metadata_tsv
)

# Replace taxonomy table
TAX <- tax_table(as.matrix(Taxonomy))
physeq <- merge_phyloseq(physeq, TAX)

# Basic sanity checks
cat("Samples:", nsamples(physeq), "\n")
cat("ASVs:", ntaxa(physeq), "\n")

# Check data format
tax_table(physeq)[1:5, ]
otu_table(physeq)[1:5, 1:5]

# Count of genera
length(get_taxa_unique(physeq, "Genus"))

# Subset to bacteria only
bac_physeq=subset_taxa(physeq, Kingdom=="Bacteria")

cat("Bacterial ASVs:", ntaxa(bac_physeq), "\n")

# Change ASV names to simpler naming convention
# Get current ASV names
asv_names <- taxa_names(bac_physeq)
# Create new names: ASV1, ASV2, ...
new_asv_names <- paste0("ASV", seq_along(asv_names))
# Create a named vector: names are old, values are new
asv_renaming <- setNames(new_asv_names, asv_names)
# Rename taxa
taxa_names(bac_physeq) <- asv_renaming
Taxonomy <- data.frame(tax_table(bac_physeq))

# Remove singletons (there shouldn't be any from dada2 but it is worth checking) 
physeq2 <- prune_taxa(taxa_sums(bac_physeq) > 1, bac_physeq)
physeq2 <- prune_samples(sample_sums(bac_physeq) > 0, bac_physeq)

# Confirm no singletons
cat("ASVs after singleton removal:", ntaxa(physeq2), "\n")
# Now revert to the bac_physeq object
bac_physeq <- physeq2
rm(physeq2)
###################################################################
# SECTION 2: Initial quality control and prevalence diagnostics
####################################################################
# NOTE: No prevalence-based filtering is applied by default in this section.
# This section is intended for exploratory QC and informed decision-making.

# Before we look at negative controls in section 3 - we do some light filtering
# Remove chloroplasts and mitochondria in one command
ps_filt <- subset_taxa(
  bac_physeq,
  (is.na(Order) | Order != "Chloroplast") &
    (is.na(Family) | Family != "Mitochondria") &
    (is.na(Genus) | Genus != "Mitochondria")
)

# This part of the code is taken from the phyloseq f1000 workflow
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html 
# Full acknowledgement to their wonderful work here

# View a table of ps_filt containing the number of features/ASVs for each phyla
# phyla with values of 1 indicate only one feature/ASV has been observed for that phyla - values here are worth noting for later
table(tax_table(ps_filt)[, "Phylum"], exclude = NULL)

# Remove any phyla characterised as NA
ps0 <- subset_taxa(ps_filt, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute the prevalence of each feature/ASV and store as a data.frame
# This command calculates prevalence = number of samples in which each ASV is present 
# (i.e., has non-zero counts).
# Resulting object is prevdf: a vector with one value per ASV.

prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Now we create a data.frame where each row is an ASV and columns include:
# Prevalence = number of samples where ASV is present
# TotalAbundance = total read count across all samples
# tax_table(ps0) = adds taxonomy data for each ASV

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

# Compute both the total and average prevalences of the features in each phyla

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# This groups by Phylum and gives:
# Mean prevalence across all ASVs in that phylum (i.e., how often ASVs 
# from this phylum appear in samples on average)
# Total prevalence = sum of all prevalence values for that phylum (i.e., total
# number of ASV appearances across all samples)
# This allows us to do prevalence filtering.

### IMPORTANT - AT THIS POINT YOU SHOULD MANUALLY INSPECT THE OUTPUT AND FILTER BASED ON THE VALUES YOU SEE.
# For example a phyla that contains average prevalence of 1 and total prevalence of 1 could reasonably be pruned 
# However, this is your dataset - with your own sample numbers and biological plausibility is important to note here too.
# You can always choose to skip this step and keep taxa to explore in more detail later in the workflow.

# If choosing to filter low prevalence phyla, this is where you would define them
# Change these to match your data if necessary, based on counts and prevalence
# The following code is commented out, as it serves as an example only
# filterPhyla = c("Phyla1", "Phyla2")
## Filter the phyloseq object to remove the specified phyla
# ps0 = subset_taxa(ps0, !Phylum %in% filterPhyla)
# ps0

# Now make a note of the number of features/ASVs and the number of genera post light filtering.
ntaxa(ps0)
length(get_taxa_unique(ps0, "Genus"))

###################################################################
# SECTION 3: Characterisation and removal of contaminants
##################################################################

# This section of the workflow directly references the decontam program and associated vignette:
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# Create a dataframe of phyloseq object in order to plot library size against sample type

df <- as.data.frame(sample_data(ps0))

# Use the unfiltered dataset to show library size

df$LibrarySize <- sample_sums(ps0)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# In this plot we are indicating the different sample-types by shape. 
# You should amend the code if you wish to plot any controls, blanks, or environmental samples included in your dataset.

ggplot(data = df, aes(x = Index, y = LibrarySize, shape = sampletype)) +
  geom_point(aes(), size = 2.5, alpha = 0.8) +
  labs(title = "Library Size against Samples",
       shape = "Sample Type") +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(face = "bold", angle = 90),
        plot.title = element_text(size = 14, face = "bold"))

# There are two options to identify contaminant taxa with the decontam package. 
# The most reliable is to use the prevalence method - based on neg control sample / blank
# Principle here is that prevalence (presence/absence across samples)
# of each feature in true samples is compared to prevalence in blanks
# to id contaminants.

# First summarise the group variable as a logical
### IMPORTANT NOTE
# It is good practice to use more than one blank/negative control. Use of a single blank, even in a small
# dataset, results in a statistically coarse assessment of prevalence. 

sample_data(ps0)$is.neg <- sample_data(ps0)[[control_column]] %in% negative_control_labels
contamdf.prev <- isContaminant(ps0,
                               method = "prevalence",
                               neg = "is.neg")
table(contamdf.prev$contaminant)
# Inspection of the table should tell you how many taxa have been flagged as potential contaminants.

head(which(contamdf.prev$contaminant))
# This code will give you the numbers of the relevant taxa. 
# You can try some more aggressive settings by setting a p value at 0.5 rather than 0.1
# this threshold will id sequences that are more prevalent in neg controls
# than in true samples, as contaminant.
contamdf.prev05 <- isContaminant(ps0,
                                 method = "prevalence",
                                 neg = "is.neg",
                                 threshold = 0.5)
table(contamdf.prev05$contaminant)
head(which(contamdf.prev05$contaminant))
# Again this code will give you a list of ASVs and you should inspect the taxonomy.
# At this point you may wish to inspect the taxonomy of these ASVs.
# I recommend inspecting, as in low-biomass samples, it is possible that taxa from samples of higher
# microbial biomass can hop across, within your dataset. I have seen this occur between a DNA mock and a blank,
# or a high depth biological sample into the blank. So it is always worth reviewing.

# Get ASV names of contaminants
contaminant_asvs_05 <- rownames(contamdf.prev05)[contamdf.prev05$contaminant]
# Subset taxonomy table to just contaminants
taxonomy_table <- Taxonomy[contaminant_asvs_05, ]
# View the result
taxonomy_table
# Now review those contaminants identified using a lower prevalence threshold
contaminant_asvs_01 <- rownames(contamdf.prev)[contamdf.prev$contaminant]
taxonomy_table_01 <- Taxonomy[contaminant_asvs_01, ]
taxonomy_table_01

# How many times were these taxa observed in all samples?
physeq.pa <- transform_sample_counts(ps0,
                                     function(abund) 1*(abund>0))
physeq.pa.neg <- prune_samples(sample_data(physeq.pa)$sampletype == "blank",
                               physeq.pa)
physeq.pa.pos <- prune_samples(sample_data(physeq.pa)$sampletype != "blank",
                               physeq.pa)
# Make a data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg),
                    contaminant=contamdf.prev05$contaminant)
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point() +
  xlab("Prevalence (Negative Controls)") +
  ylab("Prevalence (True Samples)") +
  labs(
    title = "Prevalence of ASVs in Negative Controls vs. Biological Samples",
    subtitle = "Blue = contaminants; Red = non-contaminants",
    caption = “This plot allows inspection of prevalence differences between negative controls and biological samples”,
    colour = "Contaminant")
                                     
# This plot should give you an indication as to how often the taxa identified as
# true contaminants are present within the biological samples.

## The following section gives an example of exploratory analysis, should you suspect,
# based on taxonomy, whether index hopping has occured from one of your own samples, into the blank.
# As this section is optional - the following code has been commented out
                                     
# Additional pkgs
# library(forcats); library(stringr); library(scales)

# Inputs you may need to tweak
# grp_col <- "sampletype"                 
# sample_data column with "Sample"/"Blank"
# lvl_sample <- "sample"; lvl_blank <- "blank"

# short whitelist containing taxa of interest
# wl_regex <- regex("^(TaxaName1|TaxaName2|TaxaName3|TaxaName3)$", ignore_case = TRUE)
# Or use ASV numbers
# wl_regex <- regex("^(ASV81|ASV82|ASV79|ASV44|ASV85|ASV87|ASV47|ASV45|ASV4|ASV86)$", ignore_case = TRUE)

# calculate RA at Genus level
# ps_rel <- transform_sample_counts(ps0, function(x) x/sum(x))
# gen <- tax_glom(ps_rel, taxrank = "Genus", NArm = TRUE)
# df  <- psmelt(gen) |>
#  mutate(Genus = as.character(Genus),
#         Type = ifelse(group == "blank", "blank", "sample")) |>
#  filter(str_detect(Genus, wl_regex))
# --- 1) RA for ASVs ---
# blankdf <- psmelt(ps_rel) |>
#  filter(sampletype == "blank",
#         str_detect(OTU, wl_regex))
# biodf <-   psmelt(ps_rel) |>
#  filter(sampletype != "blank",
#         str_detect(OTU, wl_regex))
                                  
# Create Summaries: mean RA, prevalence, sample enrichment
# eps <- 1e-6
# summ <- df |>
#  group_by(Genus, sampletype) |>
#  summarise(meanRA = mean(Abundance, na.rm=TRUE),
#            prev   = mean(Abundance > 0, na.rm=TRUE),
#            .groups = "drop_last") |>
#  tidyr::pivot_wider(names_from = sampletype, values_from = c(meanRA, prev), values_fill = 0) |>
#  ungroup() |>
#  transmute(
#    Genus,
#    log2_enrich = log2((get(paste0("meanRA_", lvl_sample)) + eps) /
#                         (get(paste0("meanRA_", lvl_blank )) + eps)),
#    prev_sample = get(paste0("prev_",   lvl_sample)),
#    prev_blank  = get(paste0("prev_",   lvl_blank)),
#    prev_diff   = prev_sample - prev_blank
#  )
# blank_summ <- blankdf |>
#  group_by(OTU) |>
#  summarise(meanRA = mean(Abundance, na.rm=TRUE),
#            maxRA = max(Abundance, na.rm=TRUE),
#            prev   = mean(Abundance > 0, na.rm=TRUE),
#            .groups = "drop_last")
# bio_summ <- biodf |>
#  group_by(OTU) |>
#  summarise(meanRA = mean(Abundance, na.rm=TRUE),
#            maxRA = max(Abundance, na.rm=TRUE),
#            prev   = mean(Abundance > 0, na.rm=TRUE),
#            .groups = "drop_last")
                                  
# Create a Dotplot: enrichment (y) with prevalence cue
# ggplot(summ, aes(x = fct_reorder(Genus, log2_enrich), y = log2_enrich)) +
#  geom_hline(yintercept = 0, linetype = 2) +
#  geom_point(aes(size = prev_sample, colour = log2_enrich)) +
#  coord_flip() +
#  scale_y_continuous(name = "log2(mean Relative abundance: Sample / Blank)") +
#  scale_size(range = c(2,6), name = "Prevalence in samples") +
#  scale_color_gradient2(low = "steelblue", mid = "grey60", high = "firebrick", midpoint = 0,
#                        name = "log2(mean RA Sample / mean RA Blank)") +
#  labs(x = NULL, title = "Log2 fold-change of mean relative abundance of key taxa (sample vs blank)") +
#  theme_bw(base_size = 10) +
#  theme(text = element_text(size = 10),
#        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
#        legend.background = element_blank(),
#        legend.box.background = element_rect(colour = "black"),
#        axis.text.x = element_text(face = "bold", angle = 90),
#        plot.title = element_text(size = 13, face = "bold"))

# Inspection of bio_summ may reveal ASVs that have prevalence of 0 - these you can safely remove
# Example given below.
# badtaxa <- c("ASV80", "ASV82", "ASV86", "ASV87")

# Comparison of the two _summ objects may also reveal ASVs that have a prevalence at least 10x higher in bio samples
# It is recommended you keep these taxa
# goodTaxa <- setdiff(taxa_names(ps0), badtaxa)
                                  
# ps_filt <- prune_taxa(goodTaxa, ps_filt)
# Compare the number of taxa removed
# ntaxa(ps0)
# ntaxa(ps_filt)

# Calculate RA at Genus level in the new filtered phyloseq
# ps_rel2 <- transform_sample_counts(ps_filt, function(x) x/sum(x))
# gen <- tax_glom(ps_rel2, taxrank = "Genus", NArm = TRUE)
# df  <- psmelt(gen) |>
#  mutate(Genus = as.character(Genus),
#         Type = ifelse(group == "blank", "blank", "sample")) |>
#  filter(str_detect(Genus, wl_regex))

# Summaries: mean RA, prevalence, sample enrichment
# summ <- df |>
#  group_by(Genus, sampletype) |>
#  summarise(meanRA = mean(Abundance, na.rm=TRUE),
#            prev   = mean(Abundance > 0, na.rm=TRUE),
#            .groups = "drop_last") |>
#  tidyr::pivot_wider(names_from = sampletype, values_from = c(meanRA, prev), values_fill = 0) |>
#  ungroup() |>
#  transmute(
#    Genus,
#    log2_enrich = log2((get(paste0("meanRA_", lvl_sample)) + eps) /
#                         (get(paste0("meanRA_", lvl_blank )) + eps)),
#    prev_sample = get(paste0("prev_",   lvl_sample)),
#    prev_blank  = get(paste0("prev_",   lvl_blank)),
#    prev_diff   = prev_sample - prev_blank
#  )

# Dotplot: enrichment (y) with prevalence cue
# ggplot(summ, aes(x = fct_reorder(Genus, log2_enrich), y = log2_enrich)) +
#  geom_hline(yintercept = 0, linetype = 2) +
#  geom_point(aes(size = prev_sample, colour = log2_enrich)) +
# coord_flip() +
#  scale_y_continuous(name = "log2(mean Relative abundance: Sample / Blank)") +
#  scale_size(range = c(2,6), name = "Prevalence in samples") +
#  scale_color_gradient2(low = "steelblue", mid = "grey60", high = "firebrick", midpoint = 0,
#                        name = "log2(mean RA Sample / mean RA Blank)") +
#  labs(x = NULL, title = "Log2 fold-change of mean relative abundance of key taxa (sample vs blank)") +
#  theme_bw(base_size = 10) +
#  theme(text = element_text(size = 10),
#        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
#        legend.background = element_blank(),
#        legend.box.background = element_rect(colour = "black"),
#        axis.text.x = element_text(face = "bold", angle = 90),
#        plot.title = element_text(size = 13, face = "bold"))
                                   
# summ |> arrange(desc(log2_enrich))

# Finally prune the ASVs that you have determined to be contaminants based on the above analysis                                     
# You can filter out any taxa you don't wish to remove
# using the following code example

contamdf.prev05 <- contamdf.prev05 |>
  tibble::rownames_to_column("ASV") |>
  filter(ASV != "ASVtokeep") |>
  tibble::column_to_rownames("ASV")

contam_taxa <- rownames(contamdf.prev05)[contamdf.prev05$contaminant == TRUE]
ps_decontam <- prune_taxa(!(taxa_names(ps0) %in% contam_taxa), ps0)
ntaxa(ps0)
ntaxa(ps_decontam)

# Finally, there is also the possibility to try using quantitative data.
# However, when working with host-associated microbiota - it is not possible to use any
# quantitative data based on DNA concentration (initial Qubit data) of your samples, as host DNA is likely
# to represent a large part of this number
                                     
# An alterntative is to use the post PCR input DNA concentration as provided
# by a sequencing provider, or your own work.
# NB - therefore the following data should be interpreted cautiously and not used as the sole method for
# determining contaminant taxa.

contamdf.freq <- isContaminant(ps_decontam, method="frequency", conc="quant_reading")
head(contamdf.freq)
# This calculation has returned a data.frame with several columns
# the most important being $p which contains the probability 
# that was used for classifying contaminants, and $contaminant 
# which contains TRUE/FALSE classification values with TRUE 
# indicating that the statistical evidence that the associated 
# sequence feature is a contaminant exceeds the user-settable 
# threshold. As we did not specify the threshold, the default 
# value of threshold = 0.1 was used, and $contaminant=TRUE if 
# $p < 0.1.
table(contamdf.freq$contaminant)

# Get a logical vector of contaminants
contaminant_logical <- contamdf.freq$contaminant
# Get ASV IDs (these are the rownames of contamdf.freq,
# which should match taxa_names in your phyloseq object)
contaminant_ASVs <- rownames(contamdf.freq)[contaminant_logical]
noncontaminant_ASVs <- rownames(contamdf.freq)[!contamdf.freq$contaminant]
# This will print the ASV IDs classified as contaminants
contaminant_taxonomy <- Taxonomy[contaminant_ASVs, ]

# After reviewing the taxonomy of the above, you may choose to remove them
ps_decontam
ps_decontam_freq <- prune_taxa(!contamdf.freq$contaminant, ps_decontam)
ntaxa(ps_decontam_freq)
length(get_taxa_unique(ps_decontam_freq, "Genus"))

#######################################################################
# SECTION 4: Analysis of mock samples
#######################################################################
# For this part of the analysis, you will need to have a csv file containing the expected composition of any mock communities that
# you have included in your run. Important: Your csv file should contain these fields, which will allow you to bind with your phyloseq object.
# (OTU, Sample, .data[[control_column]], Abundance, Kingdom, Phylum, Class, Order, Family, Genus, Species)
# In the example below, the Sample field is given a value of `Expected` 
# As stated earlier, you can use a combination of cell and DNA mocks.
# The code block provided allows you to visualise expected v observed relative abundance (RA) in mock samples, but you will
# need to amend the code block, or add to it, if you have used more than one type of mock (as is recommended).
# I have also included a block that focuses on analysis of logarithmically distributed mocks.

barplot_threshold <- 0.01
scatter_threshold <- 0.001
log_mock_threshold <- 1e-10
                                     
# First we create some phyloseq objects where counts are transformed into relative abundance (RA)
psrelabund <- transform_sample_counts(ps_decontam, function(x) x/sum(x))
psrelabund_unfiltered <- transform_sample_counts(bac_physeq, function(x) x/sum(x))
# 4.1 Expected vs observed composition (barplots)
# 1. Melt the phyloseq object to a long dataframe
df <- psmelt(psrelabund)

# 2. Select only your mock samples
mock_df <- df |>
  filter(.data[[control_column]] %in% mock_labels) |>
  select(OTU, Sample, .data[[control_column]], Abundance, Kingdom, Phylum, Class, Order, Family, Genus, Species)

# 3. Summarize ASV abundances across your samples (per sample)
mock_asv_table <- mock_df |>
  group_by(Sample, OTU, Genus, Species) |>
  summarise(Abundance = sum(Abundance)) |>
  arrange(Sample, desc(Abundance))

# 4. View the table
print(mock_asv_table)

write.csv(
  mock_asv_table,
  file.path(output_dir, "mock_asv_table.csv"),
  row.names = FALSE
)

# Read in input values for Mock - i.e. expected abundance
Mock_Expected_RA <- readr::read_csv(mock_expected_path)
expected_taxa_standard <- unique(Mock_Expected_RA$Genus)
print(expected_taxa_standard) # to check if all looks correct
                                                 
#Join the two dataframes
Join <- dplyr::bind_rows(
  mock_df |> mutate(Source = "Observed"),
  Mock_Expected_RA |> mutate(Source = "Expected")
)
Join <- Join |>
  mutate(Genus = ifelse(Abundance < barplot_threshold, "Other", Genus))
# Expected taxa below the plotting threshold are grouped as ‘Other’ and treated as unexpected for visualisation.
Join <- Join |>
  mutate(Expected = ifelse(Genus %in% expected_taxa_standard, "Expected", "Unexpected"))
# List your genera as they appear in the plot
genera_levels <- unique(Join$Genus)

# This produces a barplot of relative abundance of expected and observed side by side
Join |>
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "Stack", show.legend = TRUE,  width = 0.3, colour = "black") +
  scale_fill_igv() +
  labs(title = "Composition of Mock Community Standards against Expected RA",
       subtitle = "Data has been filtered to remove spurious ASVs below 0.01% relative abundance",
       x = "Sample",
       y = "Relative Abundance") +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 10, face = "bold", angle = 90),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10))
                                                 
# Analysis of your mock is another useful method to detect contaminant taxa.
# It may be worth running this code on both filtered and unfiltered phyloseq objects, to see if there is any difference.
# in the contaminant taxa that you have detected.

# 4.2 Detection thresholds from dilution mocks
# If you have run a dilution series of mock samples (also highly recommended)
# then you can use these data to determine a minimum detection threshold for expected ASVs in these samples
min_detected <- min(mock_df$Abundance[
  mock_df$Genus %in% expected_taxa_standard &
    mock_df$Abundance > 0
])
min_detected
# This should give you the min abundance which you can use to determine a filtering threshold in order to remove 
# any spurious ASVs with abundance below this level.
# For a more detailed breakdown of abundances of spurious versus expected taxa try this.
# 1. Minimum abundance of expected taxa per sample
min_expected <- Join |>
  filter(Expected == "Expected") %>%
  group_by(Sample) %>%
  dplyr::summarise(
    min_expected_abundance = min(Abundance),
    .groups = "drop"
  )
# 2. Count and maximum abundance of spurious taxa ("Other") per sample
spurious_stats <- Join %>%
  filter(Expected == "Unexpected") %>%
  filter(Abundance > 0) |>
  group_by(Sample) %>%
  dplyr::summarise(
    n_spurious_taxa = n(),
    max_spurious_abundance = max(Abundance),
    .groups = "drop"
  )
# 3. Combine both
mock_summary <- full_join(min_expected, spurious_stats, by = "Sample")
mock_summary
                                                 
# 4.3 Observed vs expected (log–log scatter)
# We can also analyse mock samples by producing a log-log observed v expected scatterplot.                                              

# Helper function to harmonise genus labels between expected and observed data
#
# This is mainly useful when a known mock taxon has only been classified
# at a higher taxonomic rank in your dataset.
#
# Example shown here:
# In some datasets, Micrococcus ASVs are only classified to Order level
# (Micrococcales). This function maps those ASVs back to "Micrococcus"
# so that expected vs observed comparisons can be made at genus level.
#
# IMPORTANT:
# You should adapt or extend this function based on the taxonomic
# behaviour observed in *your own* dataset.

canonical_genus <- function(Genus, Family = NULL, Order = NULL) {

  g <- Genus

  # Rescue Micrococcus classified only at Order level
  g <- ifelse(
    is.na(g) & !is.na(Order) & Order == "Micrococcales",
    "Micrococcus",
    g
  )

  # Clean up whitespace and missing values
  g <- ifelse(!is.na(g), stringr::str_squish(g), g)
  g[is.na(g) | g == ""] <- "Unassigned"

  return(g)
}
# Additional examples you *may* wish to include, depending on your dataset:
#
# # Harmonise Escherichia / Shigella labels
# g <- ifelse(
#   !is.na(g) & g %in% c("Escherichia", "Shigella", "Escherichia-Shigella"),
#   "Escherichia-Shigella",
#   g
# )

# OBSERVED (MC) 
obs_df <- mock_df |>
  dplyr::mutate(Canonical = canonical_genus(Genus, Family, Order)) |>
  dplyr::group_by(Sample, Canonical) |>
  dplyr::summarise(obs = sum(Abundance, na.rm = TRUE), .groups = "drop")

# EXPECTED 
exp_df <- Mock_Expected_RA |>
  dplyr::mutate(Canonical = canonical_genus(Genus, Family, Order)) |>
  dplyr::group_by(Canonical) |>
  dplyr::summarise(exp = sum(Abundance, na.rm = TRUE), .groups = "drop")

# JOIN BY TAXON ONLY
comp <- full_join(exp_df, obs_df, by = "Canonical") |>
  replace_na(list(exp = 0, obs = 0))

# define the mock taxa once from expected
mock_taxa <- exp_df$Canonical

# Optional: filter to keep all mock taxa + non-mock >= 0.1% (0.001 threshold)
# Threshold used only to collapse low-abundance taxa for visual clarity
comp_filt <- comp |>
  filter(Canonical %in% mock_taxa | pmax(exp, obs) >= scatter_threshold)

# flags for plotting
comp_filt <- comp_filt |>
  dplyr::mutate(
    absent_in_obs  = exp > 0 & obs == 0,
    unexpected_tax = !(Canonical %in% mock_taxa) & obs > 0
  )
# SCATTERPlot: observed vs expected (single panel)
eps <- 1e-4
# make a flag for mock taxa
comp3 <- comp_filt |>
  mutate(is_mock = exp > 0)

# per-facet metrics on the mock taxa only (exp > 0) ---
eps_stat <- 1e-6
stats_df <- comp3 |>
  dplyr::filter(is_mock) |>
  mutate(exp_adj = pmax(exp, eps_stat),
         obs_adj = pmax(obs, eps_stat),
         lx = log10(exp_adj), ly = log10(obs_adj),
         log_bias = ly - lx,
         abs_log_err = abs(log_bias)) |>
  group_by(Sample) |>
  group_modify(~{
    m <- lm(ly ~ lx, data = .x)
    tibble(
      n         = nrow(.x),
      slope     = as.numeric(coef(m)[2]),
      r2        = summary(m)$r.squared,
      med_bias  = median(.x$log_bias),
      mae       = median(.x$abs_log_err)
    )
  }) |>
  ungroup() |>
  # round and build facet-strip labels with sample info
  mutate(
    slope    = round(slope, 2),
    r2       = round(r2, 2),
    med_bias = round(med_bias, 2),
    mae      = round(mae, 2)
  )

# plot
p1 <- ggplot(comp3, aes(pmax(exp, eps), pmax(obs, eps))) +
  # 1:1 and tolerance bands (±2× dotted, ±5× dotted lighter)
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5) +
  geom_abline(slope = 2, intercept = 0, linetype = "dotted", alpha = 0.6) +
  geom_abline(slope = 0.5, intercept = 0, linetype = "dotted", alpha = 0.6) +
  geom_abline(slope = 5, intercept = 0, linetype = "dotted", alpha = 0.35) +
  geom_abline(slope = 0.2, intercept = 0, linetype = "dotted", alpha = 0.35) +
  # unexpected taxa (non-mock): faint solid dots
  geom_point(data = subset(comp3, !is_mock),
             size = 1.3, shape = 16, alpha = 0.5) +
  # expected taxa (mock): hollow circles
  geom_point(data = subset(comp3, is_mock),
             size = 2.8, shape = 21, fill = "white", stroke = 0.9) +
  # labels only for expected taxa
  ggrepel::geom_text_repel(data = subset(comp3, is_mock),
                           aes(label = Canonical),
                           size = 2, max.overlaps = Inf) +
  facet_wrap(~ Sample, nrow = 1) +
  scale_x_log10(limits = c(1e-4, 1), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_log10(limits = c(1e-4, 1), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  labs(
    title = "Scatterplot of observed vs expected taxa in mock samples (genus level)",
    subtitle = "Dashed = 1:1; dotted = ±2× and ±5× bands. Hollow = mock taxa; faint = unexpected (≥0.1% RA).",
    x = "log10 Expected Relative Abundance",
    y = "log10 Observed Relative Abundance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(face = "bold", angle = 90, vjust = 0.5),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "none",
    plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")
  ) + annotate("text", x = 7e-4, y = 0.030, label = "±2×", size = 3) +
  annotate("text", x = 7e-4, y = 0.240, label = "±5×", size = 3)

p1
print(stats_df)
# This last command allows you to inspect the regression metrics you have calculated
# slope, R^2, median log bias, MAE.

# 4.4 Logarithmic mock standards 
# This last code block is used for any mocks with a logarithmic distribution of taxa.
# The code assumes your mock is named LogMC
# You should create a csv file of the expected abundances as before.
# Read in input values for Mock - i.e. expected abundance
Log_Mock_RA <- readr::read_csv(log_mock_expected_path)
expected_taxa_log <- unique(Log_Mock_RA$Genus)
print(expected_taxa_log) # check if all good

# Reset the data
log_mock_df <- psmelt(psrelabund)

# Filter for the log standards
log_mock_df <- log_mock_df |>
  filter(Sample == log_mock_sample) |>
  select(OTU, Sample, .data[[control_column]], Abundance, Kingdom, Phylum, Class, Order, Family, Genus, Species)
Join <- rbind(log_mock_df, Log_Mock_RA)
Join <- Join |>
  mutate(Genus = ifelse(Abundance < log_mock_threshold, "Other", Genus))
Join |>
  ggplot(aes(x = Sample, y = reorder(Genus, +Abundance), size = Abundance, color = Genus)) +
  geom_point(alpha = 0.8) +
  scale_colour_igv() +
  scale_size_continuous(trans = "log10") +
  facet_wrap(~.data[[control_column]], scales = "free_x") +
  labs(title = "Abundance of Genera in a Logarithmic DNA Mock Standard",
       caption = "Data has been filtered to remove spurious ASVs below 0.000001% relative abundance",
       y = "Genus", x = "Sample", size = "Abundance (log10)") +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(color = "none")

# As with dilution series of mock samples, the log sample can be used to determine the minimum detection threshold level
# this can be used as a defensible filtering level, at the ASVs at lower abundances are likely to be spurious.
min_detected_dna <- min(log_mock_df$Abundance[
  log_mock_df$Genus %in% expected_taxa_log &
    log_mock_df$Abundance > 0
])
min_detected_dna
#######################################
# SECTION 5: Final Prevalence Filtering
########################################

# If you have included mocks in your run then you can use the analysis above to determine the minimum detection threshold.
# For the purposes of this workflow I have specified a threshold based on the minimum detection thresholds calcyulated from
# one of my own datasets, but you should amend this figure to reflect your own dataset.
mock_abund_threshold <- 0.00014
# Use this threshold to create a filtered dataframe of all samples and filtered phyloseq objects.
filtered_df <- df[df$Abundance >= mock_abund_threshold, ]
# Identify taxa to keep
keep_taxa <- taxa_names(psrelabund)[taxa_sums(psrelabund) >= mock_abund_threshold * nsamples(psrelabund)]
ps_filtered <- prune_taxa(keep_taxa, ps_decontam)
ps_rel_filtered <- transform_sample_counts(ps_filtered, function(x) x / sum(x))
# Compare number of ASVs before and after filtering
ntaxa(bac_physeq)
ntaxa(ps_decontam)
ntaxa(ps_filtered)
# Note number of genera
length(get_taxa_unique(ps_filtered, "Genus"))
# Let's visualise the phyla.
# Result: a vector with one value per ASV.
prevdf2 = apply(X = otu_table(ps_filtered),
               MARGIN = ifelse(taxa_are_rows(ps_filtered), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Now you create a data.frame where each row is an ASV and columns include:
# Prevalence = number of samples where ASV is present
# TotalAbundance = total read count across all samples
# Add taxonomy and total read counts to this data.frame
prevdf2 = data.frame(Prevalence = prevdf2,
                    TotalAbundance = taxa_sums(ps_filtered),
                    tax_table(ps_filtered))
# Compute the total and average prevalences of the features in each phyla
plyr::ddply(prevdf2, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Subset to the remaining phyla
prevdf3 = subset(prevdf2, Phylum %in% get_taxa_unique(ps_filtered, "Phylum"))
ggplot(prevdf3, aes(TotalAbundance, Prevalence / nsamples(ps_filtered),color=Phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) +
  labs(title = "Prevalence versus total counts of ASVs",
       subtitle = "Dataset has been filtered to remove contaminants and low prevalence ASVs") +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        axis.text.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "none")
# Combine the neg and mock control labels
control_labels <- c(mock_labels, negative_control_labels)
# Use this to subset
final_physeq <- subset_samples(
  ps_filtered,
  !(get(control_column) %in% control_labels)
)
# Sanity check
nsamples(final_physeq)
sample_names(final_physeq)
# Repeat for the RA object
final_rel_physeq <- subset_samples(
  ps_rel_filtered,
  !(get(control_column) %in% control_labels)
)
filtered_df <- filtered_df |>
  dplyr::filter(!(.data[[control_column]] %in% control_labels))

# Inspect the counts in the final phyloseq object
nb_samples  <- nsamples(final_physeq)
nb_samples
nb_features <- ntaxa(final_physeq)
nb_features
length(get_taxa_unique(genus_physeq, "Genus"))
# This should tell you whether pruning contaminants, filtering, removing mocks and the blank has worked

# Following filtering you should consider zero ASV counts - i.e. a measure of sparsity in the data-set
sum(as(otu_table(final_physeq), "matrix") == 0)
sum(as(otu_table(final_physeq), "matrix") == 0) / (nb_features * nb_samples) * 100
# This shows number of zero ASV counts (ASVs that do not appear in a sample)
# And the percentage represents the proportion of ASVs in the dataset are zero count ASVs                                           

# A zero may mean that this feature is not present in a sample,
# Or it may mean that sequencing depth was inadequate

# Create a matrix of the otu table
mat <- t(otu_table(final_physeq))
class(mat) <- "matrix"
class(mat)
#> [1] "matrix" "array"

# What is the minimum total count over all samples per OTU in this data set?
min(colSums(mat))
non_zero <- 0*1:nb_features

for (i in 1:nb_features){
  non_zero[i]<-sum(final_physeq@otu_table[i, ] != 0)
}

plot(non_zero, xlab = "Unique ASVs", ylab = "Frequency of Samples (n)", 
     main = "Frequency of Non-Zero ASV occurrences in dataset, by sample.", las = 1)

# Generate rarefaction curves with ggrare
library(ranacapa)
rarecount <- ggrare(final_physeq, step = 100, se = FALSE) +
  # geom_vline(xintercept = 2000, linetype = "dashed", color = "orange") + # use this if want to display a potential depth-threshold
  labs(
    title = "Alpha Rarefaction Curves",
    subtitle = "Rarefied at incremental steps; vertical line = 2000 read cutoff",
    x = "Sequencing Depth",
    y = "Observed ASVs"
  ) +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 10),
    plot.margin = ggplot2::margin(1, 1, 0.5, 1, "cm"),
    plot.subtitle = element_text(size = 7),
    axis.text = element_text(size = 7, angle = 90),
    legend.text = element_text(size = 6)
  )
rarecount
# This plot should illustrate whether you have achieved suitable sequencing depth in the samples.
# The count at which the samples curves plateau is considered to be a suitable depth - and could be considered the lowest threshold for depth.
# Make a note, replot with the geom_vline command amended to include your threshold. This number should inform the next plot.

# Now to plot sample counts
sum_seq <- rowSums(mat)
plot(sum_seq, ylim=c(0,60000), main=c("Number of counts per sample"), 
     xlab=c("Samples"), ylab = c("Read count(n)"))
sum_seq
min(sum_seq)
max(sum_seq)
# This is a quick visualisation to check whether any of the samples have drastically been reduced in read size 
# following filtering and decontamination.
# This is especially important in host-associated low-microbial-biomass datasets which mayu be low depth to begin with.

# Data needs to be normalised in future analysis (log for beta, relative abundance for bar plots)
# Let's make some more visually appealing plots
sample_sums_df <- data.frame(SampleID = sample_names(final_physeq),
                             ReadCount = sample_sums(final_physeq))
histcount <- ggplot(sample_sums_df, aes(x = ReadCount)) +
  geom_histogram(aes(),bins = 30, color = "black") +
  # geom_vline(xintercept = 2000, linetype = "dashed", color = "black") + # use this to indicate a potential cut-off for sequencing depth
  # annotate("text", x = 2100, y = 11.5, label = "Rarefaction threshold (2000)", hjust = 0, size = 2) +
  labs(title = "Distribution of Sample Read Counts",
       subtitle = "X intercept represents sample depth threshold",
       x = "Total ASV Counts per Sample", y = "Number of Samples") +
  scale_fill_npg() +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 7),
        plot.margin = ggplot2::margin(0.5, 1, 1, 1, "cm"),
        axis.text.x = element_text(size = 6, face = "bold", angle = 90),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7),
        legend.text = element_text(size = 6))
histcount
# Overall this combination of plots should tell you whether you have any samples that should be pruned from the dataset owning to low depth
# And your histogram of sample counts will give you an idea of how many samples you might loose based on the depth pruning.
# Once you have made your decision proceed with the optional pruning command:

# Remove low-depth samples up front
# low_depth_samples <- sample_names(final_physeq)[sample_sums(genus_physeq) < 2000]
# final_physeq_no_outliers <- subset_samples(final_physeq, !(sample_names(final_physeq) %in% low_depth_samples))
# Metadata_no_outliers <- Metadata %>% filter(!(sampleid %in% low_depth_samples))

# If you have pruned wisely, you will find that there is little difference in the following counts prior to pruning.
# ntaxa(final_physeq_no_outliers)
# length(get_taxa_unique(final_physeq_no_outliers, "Genus"))

                                                 
###################################################################
#
# Section 5: Alpha diversity metrics
#
###################################################################

# Alpha and Beta Diversity
# Using Vegan to caluclate alpha diversity
# Remove low-depth samples up front
low_depth_samples <- sample_names(genus_physeq)[sample_sums(genus_physeq) < 2000]
# This will remove 10 semen samples. Not ideal. Could set bar higher in terms of depth
# but there will be a trade of with a much greater number of samples being pruned.
genus_physeq_no_outliers <- subset_samples(genus_physeq, !(sample_names(genus_physeq) %in% low_depth_samples))
Metadata_no_outliers <- Metadata %>% filter(!(sampleid %in% low_depth_samples))
# Subset to semen samples only
ntaxa(genus_physeq_no_outliers)
# Still 585
length(get_taxa_unique(genus_physeq_no_outliers, "Genus"))
# Still 194
semen_physeq <- subset_samples(genus_physeq_no_outliers, sampletype == "semen")
Metadata_semen <- Metadata_no_outliers %>% filter(sampletype == "semen")
# Add a metadata column for semen concentration -binned to make it categorical
sample_data(semen_physeq)$concentration_binned <- cut(sample_data(semen_physeq)$Mojo_concentration,
                                                      breaks = c(0, 90, 250))
Metadata_semen$concentration_binned <- sample_data(semen_physeq)$concentration_binned
# Get total read counts per sample (after host removal)
read_depth <- sample_sums(semen_physeq)
read_depth
# Add to sample_data
sample_data(semen_physeq)$read_depth <- read_depth
sample_data(semen_physeq)$depth_binned <- cut(sample_data(semen_physeq)$read_depth,
                                                      breaks = c(2000, 5000, 10000, 60000))
Metadata_semen$read_depth <- sample_data(semen_physeq)$read_depth
Metadata_semen$depth_binned <- sample_data(semen_physeq)$depth_binned
# Repeat for whole dataset
# Get total read counts per sample (after host removal)
read_depth <- sample_sums(genus_physeq_no_outliers)
read_depth
# Add to sample_data
sample_data(genus_physeq_no_outliers)$read_depth <- read_depth
# Calculate alpha diversity metrics
otu <- t(otu_table(semen_physeq))
S <- specnumber(otu)
shannon <- diversity(otu, index = "shannon")
pielou <- shannon / log(S)
richness <- t(estimateR(otu))  # S.obs, Chao1, etc.
# Combine into a single data frame
alpha_df <- Metadata_semen %>%
  mutate(
    S.obs = richness[,"S.obs"],
    Chao1 = richness[,"S.chao1"],
    Shannon = shannon,
    Pielou = pielou
  )
# Check correlation with mapped_boar
# Correlations of alpha diversity with other technical covariates
ggplot(alpha_df, aes(x = mapped_boar, y = Chao1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Chao1 diversity vs Mapped to Boar Genome",
       x = "Reads mapped to boar genome (%)",
       y = "Chao1 Index")
ggsave("Chao1_boar_mapping.png", scale = 2)
cor.test(alpha_df$Chao1, alpha_df$mapped_boar)
# Pearson's product-moment correlation
# data:  alpha_df$Chao1 and alpha_df$mapped_boar
# t = -7.839, df = 28, p-value = 1.541e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.9156819 -0.6680569
# sample estimates:
#       cor 
# -0.8288405 
cor.test(alpha_df$S.obs, alpha_df$mapped_boar)
# Pearson's product-moment correlation
# data:  alpha_df$S.obs and alpha_df$mapped_boar
# t = -7.839, df = 28, p-value = 1.541e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.9156819 -0.6680569
# sample estimates:
#       cor 
# -0.8288405
# Mapped to boar is a proxy for read depth - so use one or the other in the analysis following
# Repeat above for read depth
# Correlations of alpha diversity with other technical covariates
ggplot(alpha_df, aes(x = read_depth, y = Chao1)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  labs(title = "Relationship Between Sequencing Depth and Alpha Diversity",
       subtitle = "Significant correlation between sample depth and Chao1 richness",
       x = "Sample Read Count",
       y = "Chao1 Index") +
  stat_cor(method = "pearson", label.x = 3000, label.y = 120) +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10))
ggsave("Chao1_depth.pdf", width = 8, height = 6, units = "in")
cor.test(alpha_df$Chao1, alpha_df$read_depth)
# Pearson's product-moment correlation
# data:  alpha_df$Chao1 and alpha_df$read_depth
# t = 6.3612, df = 28, p-value = 6.964e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.5650029 0.8841662
# sample estimates:
#       cor 
# 0.7687838  

# NB Chao1 and S.obs/richness are reporting the same results in this dataset because it 
# contains no singleton observations - this is a result of the dada2 denoising workflow.

# mapped_boar/read_depth and chao1 / S.obs extremely highly correlated.
P1 <- ggplot(alpha_df, aes(x = Site, y = S.obs)) +
  geom_boxplot(aes(fill = Site)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 120) +
  labs(title = 'Observed Richness (S.obs)', x = '', y = '', tag = "A") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

P2 <- ggplot(alpha_df, aes(x = Site, y = Chao1)) +
  geom_boxplot(aes(fill = Site)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 120) +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

P3 <- ggplot(alpha_df, aes(x = Site, y = Pielou$shannon)) +
  geom_boxplot(aes(fill = Site)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.7) +
  labs(title= 'Evenness (Pielou)', x= ' ', y= '', tag = "C") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

P4 <- ggplot(alpha_df, aes(x = Site, y = Shannon$shannon)) +
  geom_boxplot(aes(fill = Site)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 2.0) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

(P1 | P2) / (P3 | P4) +
  plot_annotation(
    title = "Alpha Diversity Metrics by Site",
    subtitle = "Richness, Chao1, Evenness (Pielou), and Shannon indices with Wilcoxon p-values",
    caption = "Note: S.obs and Chao1 are highly correlated with read depth and should be interpreted with caution.",
    theme = theme(text = element_text(size = 10),
                  plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 10),
                  plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))
  )
ggsave("Genus_Alpha_diversity_site_no_outliers.pdf", width = 8, height = 6, units = "in")

P1 <- ggplot(alpha_df, aes(x = Line, y = S.obs)) +
  geom_boxplot(aes(fill = Line)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 120) +
  labs(title = 'Observed Richness (S.obs)', x = '', y = '', tag = "A") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7))

P2 <- ggplot(alpha_df, aes(x = Line, y = Chao1)) +
  geom_boxplot(aes(fill = Line)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 120) +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7))

P3 <- ggplot(alpha_df, aes(x = Line, y = Pielou$shannon)) +
  geom_boxplot(aes(fill = Line)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.7) +
  labs(title= 'Evenness (Pielou)', x= ' ', y= '', tag = "C") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7))

P4 <- ggplot(alpha_df, aes(x = Line, y = Shannon$shannon)) +
  geom_boxplot(aes(fill = Line)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 2.0) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7))

(P1 | P2) / (P3 | P4) +
  plot_annotation(
    title = "Alpha Diversity Metrics by Boar Line",
    subtitle = "Richness, Chao1, Evenness (Pielou), and Shannon indices with Kruskal p-values",
    caption = "Note: S.obs and Chao1 are highly correlated with read depth and should be interpreted with caution.",
    theme = theme(text = element_text(size = 10),
                  plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 10),
                  plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")))
ggsave("Genus_Alpha_diversity_line_no_outliers.pdf", width = 8, height = 6)


P1 <- ggplot(alpha_df, aes(x = concentration_binned, y = S.obs)) +
  geom_boxplot(aes(fill = concentration_binned)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 120) +
  labs(title = 'Observed Richness (S.obs)', x = '', y = '', tag = "A") +
  theme(legend.position = "none")

P2 <- ggplot(alpha_df, aes(x = concentration_binned, y = Chao1)) +
  geom_boxplot(aes(fill = concentration_binned)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 120) +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  theme(legend.position = "none")

P3 <- ggplot(alpha_df, aes(x = concentration_binned, y = Pielou$shannon)) +
  geom_boxplot(aes(fill = concentration_binned)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.7) +
  labs(title= 'Evenness (Pielou)', x= ' ', y= '', tag = "C") +
  theme(legend.position = "none")

P4 <- ggplot(alpha_df, aes(x = concentration_binned, y = Shannon$shannon)) +
  geom_boxplot(aes(fill = concentration_binned)) +
  geom_point(size = 1, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 2.0) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  theme(legend.position = "none")

(P1 | P2) / (P3 | P4) +
  plot_annotation(
    title = "Alpha Diversity Metrics by Semen Concentration (Semen Samples Only, Outliers Removed)",
    subtitle = "Richness, Chao1, Evenness (Pielou), and Shannon indices with Wilcoxon p-values",
    caption = "Note: S.obs and Chao1 are highly correlated with mapping of reads to boar genome and should be interpreted with caution.",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12)))
ggsave("Genus_Alpha_diversity_concentration_no_outliers.png", scale = 3)

# Could potentially look at the skin and glove samples - but there are so few
# of them that statistical tests won't be robust - leave for now.

######################################################################
#
# Section 6: Beta diversity
#
######################################################################
# In this section we will be exploring beta-diversity using the datasets with outliers removed.
# Aim here is to explore beta diversity across the different niches (i.e. semen v gloves,
# boar swabs and extender), and the beta diversity within the semen samples based on the 
# semen metadata we have: Site, Boar Line, Semen concentration, Motility, Viability, plus
# checking for correlation with technical variables too. We have seen in alpha diversity
# that the percentage of reads that originally mapped to the boar genome (i.e. a proxy for
# read depth) was highly correlated with richness/chao1.

# Another priority here is to explore beta diversity using different distance indices.
# Particularly interested in those that incorporate phylogenetic information such as
# unifrac, but should also cover the other common and useful distances such as Bray Curtis
# for completeness.

# Let's look at some of the metadata variables that have continuous distribution
qplot(sample_data(semen_physeq)$Mojo_concentration, geom = "histogram") + xlab("Concentration")
# This shows a skewed guassian distribution - with peak around 50x10^6 cells/ml. The range
# upwards to 250 become sparser here.
qplot(sample_data(semen_physeq)$Total_Motility, geom = "histogram") + xlab("Motility")
# Skewed to the higher end - most samples between 80-95% showed large amount of homogeneity.
# There are only 4 samples lower than this, 2 around 70%, 1 at ~60 and the other 55%.
qplot(sample_data(semen_physeq)$Alive, geom = "histogram") + xlab("Viability")
# Viability shows peak at 70% and remainder of samples are mostly above this, up to 90%.
# There are a couple of outliers at lower measures - 45 and 60%.
qplot(sample_data(semen_physeq)$Acrosome_damaged, geom = "histogram") + xlab("Acrosome damage")
# Low counts here, as expected. Poisson distribution (I think - whatever the name is for
# A skewed distribution). Peak around 4-5%, then downward trend to 25%.
qplot(sample_data(semen_physeq)$DFI, geom = "histogram") + xlab("Acrosome damage")
# Again, low values - DFI stands for DNA fragmentation index. Most peak beneath 1%.
# highest value is 6% which is clear outlier. Main batch of values end at 2%.

# Overall semen metadata is largely homogenous with a few outliers, as is normal with 
# breeding boar, as these all have high quality semen. Concentration is probably the
# most interesting metric here.

qplot(sample_data(semen_physeq)$quant_reading, geom = "histogram") + xlab("DNA concentration")
# unremarkable. A large block of values between 8-22. One outlier at 35.
qplot(sample_data(semen_physeq)$mapped_boar, geom = "histogram") + xlab("Read mapped to boar genome %")
# A large range here, from 0 to 80%. Distribution fairly even.

# We know that I have a wide variety of read depths, with many samples at the low end.
# Let's look at the logged counts per sample.
# Remember that these physeqs have outliers removed, as in previous section.
qplot(log10(rowSums(otu_table(semen_physeq)))) +
  xlab("Logged counts-per-sample")
# This distribution suggests that log transformation is ok for initial quick
# look at data normalisation. The distribution is now much more gaussian.
qplot(log10(rowSums(otu_table(genus_physeq_no_outliers)))) +
  xlab("Logged counts-per-sample")
# Again, gaussian. Suggests this transformation is ok for beta diversity metrics that
# can be visualised with PCoA. Consider CLR for other measures.

# The total_motility from genus, is very homogenous and doesn't suggest there
# will be much use in making this a categorical variable.
# The mojo motility may be better, but what does it mean really? As semen
# analysis suggested, may just be a factor of the semen concentration
# semen concentration may be better to work with as a new factor.

sslog <- transform_sample_counts(semen_physeq, function(x) log(1 + x))
genuslog <- transform_sample_counts(genus_physeq_no_outliers, function(x) log(1 + x))

# Explore relatedness of all niches
ggord <- ordinate(genuslog, method = "PCoA", distance = "wunifrac")
ggevals <- ggord$values$Eigenvalues
ggord_df <- plot_ordination(genuslog, ggord, color = "sampletype", justDF = TRUE)
ggord_df$label <- ifelse(ggord_df$Axis.1 < -0.08, rownames(ggord_df), NA) # label only outliers
pc1_var <- round(100 * ggevals[1] / sum(ggevals), 1)
pc2_var <- round(100 * ggevals[2] / sum(ggevals), 1)
wuniall <- ggplot(ggord_df, aes(x = Axis.1, y = Axis.2, color = sampletype)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),
       title = "Weighted Unifrac PCoA by Sample Type",
       subtitle = "Semen and boar swabs are closely related, with two notable semen outliers",
       color = "Sample Type") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, na.rm = TRUE) +
  scale_color_discrete(label = c("Boar Swab", "Extender", "Glove swab", "Boar Semen")) +
  coord_fixed(sqrt(ggevals[2] / ggevals[1])) +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm"),
        plot.subtitle = element_text(size = 11),
        axis.text = element_text(size = 7, angle = 90),
        legend.text = element_text(size = 10))
wuniall
ggsave("wunifrac_all.png", width = 8, height = 6, dpi = 300)
ggsave("wunifrac_all.pdf", width = 8, height = 6)
# Here semen sample 0839K is a clear outlier.
# Weighted unifrac shows the most clear distinction between the sample types.
# The boar swab samples are most closely related to the semen samples.

# Performing a PERMANOVA.
# Checking model with all technical variables included - no significance found except for sampletype
# Removing these variables from final model to avoid overfitting.
adonis2(distance(genuslog, method = "wunifrac") ~ sampletype,
        data = as(sample_data(genuslog), "data.frame"),by = "margin")
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999

# adonis2(formula = distance(genuslog, method = "wunifrac") ~ sampletype, data = as(sample_data(genuslog), "data.frame"), by = "margin")
# Df SumOfSqs      R2      F Pr(>F)   
# sampletype  3  0.12315 0.24464 3.4547  0.018 *
# Residual    32  0.38025 0.75536                  
# Total       35  0.50304 1.00000                 
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Warning message:
#  In matrix(tree$edge[order(tree$edge[, 1]), ][, 2], byrow = TRUE,  :
#              data length [1155] is not a sub-multiple or multiple of the number of rows [578]

# Explore semen samples now - first by site
ord <- ordinate(sslog, method = "PCoA", distance = "wunifrac")
evals <- ord$values$Eigenvalues
ord_df <- plot_ordination(sslog, ord, color = "Site", justDF = TRUE)
ord_df$label <- ifelse(ord_df$Axis.1 > 0.1, rownames(ord_df), NA) # label only outliers
pc1_var <- round(100 * evals[1] / sum(evals), 1)
pc2_var <- round(100 * evals[2] / sum(evals), 1)
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),
       title = "Weighted Unifrac PCoA of Semen by Boar Stud",
       subtitle = "Willingham samples appear more closely related",
       caption = "Data has been normalised through log transformation and low depth samples removed",
       color = "Sample Type") +
  stat_ellipse() +
  ggrepel::geom_text_repel(aes(label = label), size = 3, na.rm = TRUE) +
  scale_color_observable() +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 10),
        plot.margin = ggplot2::margin(1, 1, 0.5, 1, "cm"),
        plot.subtitle = element_text(size = 7),
        axis.text = element_text(size = 7, angle = 90),
        legend.text = element_text(size = 6))
ggsave("wunifrac_site.png", scale = 2)
# Here semen sample 0839K is a clear outlier again, 0313K less so.
wuni_stud <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, colour=Site, size = read_depth)) +
  geom_point(alpha = 0.9) +
  labs(x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),
       title = "Weighted Unifrac PCoA of Semen by Boar Stud",
       subtitle = "Willingham samples appear more closely related, and higher depth samples are tightly clustered",
       caption = "Data has been normalised through log transformation and low depth samples removed",
       color = "Site",
       size = "Read Depth") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, na.rm = TRUE) +
  scale_color_observable() +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm"),
        plot.subtitle = element_text(size = 11),
        axis.text = element_text(size = 10, angle = 90),
        legend.text = element_text(size = 10))
wuni_stud
ggsave("wunifrac_stud_depth.pdf", width = 8, height = 6)
combo_wuni <- wuniall + wuni_stud+ plot_layout(nrow = 2)
combo_wuni
ggsave("combo_wuni.pdf", plot = combo_wuni, width = 9, height = 7, units = "in")
# Check stats
adonis2(distance(sslog, method = "wunifrac") ~ Site + read_depth,
        data = as(sample_data(sslog), "data.frame"),by = "margin")

# Now by boar line
ord <- ordinate(sslog, method = "PCoA", distance = "wunifrac")
evals <- ord$values$Eigenvalues
ord_df <- plot_ordination(sslog, ord, color = "Line", justDF = TRUE)
ord_df$label <- ifelse(ord_df$Axis.1 > 0.1, rownames(ord_df), NA) # label only outliers
pc1_var <- round(100 * evals[1] / sum(evals), 1)
pc2_var <- round(100 * evals[2] / sum(evals), 1)
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Line)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),
       title = "Weighted Unifrac PCoA of Semen by Boar Line",
       subtitle = "LW and LR samples appear most closely related, however there is substantial overlap between all groups",
       caption = "Data has been normalised through log transformation and low depth samples removed",
       color = "Boar Line") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, na.rm = TRUE) +
  scale_color_observable() +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  theme_bw()
ggsave("wunifrac_line.png", scale = 2)
# Appears to be some clustering but difficult to make out.

# Exploring the ordination as a convex hull plot.
# Assuming ord_df is your ordination data frame and 'Line' is your group variable
find_hull <- function(df) df[chull(df$Axis.1, df$Axis.2), ]
hulls <- ord_df %>% group_by(Line) %>% do(find_hull(.))
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Line)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_polygon(data = hulls, aes(fill = Line, group = Line), alpha = 0.2, color = NA) +
  labs(title = "PCoA with Convex Hulls by Boar Line",
       x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),) +
  theme_minimal()
# Centroids and segments
centroids <- ord_df %>%
  group_by(Line) %>%
  summarize(centroid1 = mean(Axis.1), centroid2 = mean(Axis.2))
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Line)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_point(data = centroids, aes(x = centroid1, y = centroid2), size = 5, shape = 4, color = "black") +
  geom_segment(data = merge(ord_df, centroids, by = "Line"),
               aes(xend = centroid1, yend = centroid2), alpha = 0.3) +
  labs(title = "PCoA with Group Centroids by Boar Line",
       x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),) +
  theme_minimal()
# Density contours
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Line, fill = Line)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_density_2d(aes(fill = Line), geom = "polygon", alpha = 0.2, color = NA) +
  labs(title = "PCoA with Density Contours by Boar Line") +
  theme_minimal()
# Faceting
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Line)) +
  geom_point(size = 3, alpha = 0.7) +
  facet_wrap(~ Line) +
  labs(title = "PCoA Faceted by Boar Line") +
  theme_minimal()
# I think the base PcoA is best for now, although the facet does show just how much
# overlap there is between the boar lines.

# Now let's try semen characterisitics
ord_df <- plot_ordination(sslog, ord, shape = "Site", justDF = TRUE)
ord_df$label <- ifelse(ord_df$Axis.1 > 0.1, rownames(ord_df), NA) # label only outliers
pc1_var <- round(100 * evals[1] / sum(evals), 1)
pc2_var <- round(100 * evals[2] / sum(evals), 1)
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, shape = Site, color = Mojo_concentration)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),
       title = "Weighted Unifrac PCoA of Semen by Site and Concentration",
       subtitle = "Highly concentrated samples are exclusively from Willingham",
       caption = "Data has been normalised through log transformation and low depth samples removed",
       color = "Semen Concentration (10^6)") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, na.rm = TRUE) +
  scale_color_bs5("purple") +
  stat_ellipse() +
  coord_fixed(sqrt(evals[2] / evals[1]))
ggsave("wunifrac_site_conc.png", scale = 2)
# important question to ask Craig here - how is extender added. Does it bear any relationship
# to the concentration of the original samples at all? If not it is pointless presenting
# this issue.

ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = read_depth, shape = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),
       title = "Weighted Unifrac PCoA of Semen by Site and Read depth",
       subtitle = "Higher depth samples appear tightly clustered.",
       caption = "Data has been normalised through log transformation and low depth samples removed",
       color = "Read depth") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, na.rm = TRUE) +
  scale_color_viridis_c() + 
  coord_fixed(sqrt(evals[2] / evals[1])) +
  theme_bw()
ggsave("wunifrac_site_depth.png", scale = 2)
# Now that low depth samples are removed - read depth still appears significant
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = read_depth, shape = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(x = paste0("PCoA 1 (", pc1_var, "%)"),
       y = paste0("PCoA 2 (", pc2_var, "%)"),
       title = "Weighted Unifrac PCoA of Semen by Site and Read depth",
       subtitle = "Despite removal of the lowest depth samples, a relationship between depth and distance remains apparent.",
       caption = "Data has been normalised through log transformation and low depth samples removed",
       color = "Read depth") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, na.rm = TRUE) +
  scale_color_bs5("purple") + 
  coord_fixed(sqrt(evals[2] / evals[1]))
ggsave("wunifrac_site_depth2.png", scale = 2)
# Code in case want to look at binned continuous variables.
sample_data(sslog)$concentration_binned <- cut(sample_data(sslog)$Mojo_concentration,
                                     breaks = c(0, 90, 250))
sample_data(sslog)$mapped_binned <- cut(sample_data(sslog)$mapped_boar,
                                               breaks = c(0, 25, 50, 100))

# P is no longer significant
wuni_dist <- phyloseq::distance(sslog, method = "wunifrac")
beta_disp <- betadisper(wuni_dist, sample_data(sslog)$concentration_binned)
anova(beta_disp)
# The beta-dispersion test is a way of testing whether differences in spread
# are responsible for the permanova result.

# This shows that Boston only contains low conc samples
# Plan - subset to Willingham
sslog_willingham <- subset_samples(sslog, Site == "Willingham")
# THINK ABOUT EXPLORING IN MORE DETAIL BASED ON WHAT CRAIG SAYS

# Alternative methods CLR
library(microbiome)
clr_phy <- microbiome::transform(sslog, "clr")
# calculate Euclidean distances
clr_dist <- dist(t(otu_table(clr_phy)))
# Ordinate and plot
ordination_clr <- ordinate(clr_phy, method = "PCoA", distance = clr_dist)
ord_df <- plot_ordination(clr_phy, ordination_clr, justDF = TRUE)

# Example plot
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = read_depth)) +
  geom_point(size = 3) +
  scale_colour_viridis_c() +
  labs(title = "PCoA of CLR-transformed Data (Aitchison Distance), coloured for read depth of samples",
       subtitle = "There is a correlation between sample read depth and Axis 1",
       caption = "Pearson's product-moment correlation, p = 1.621e-06",
       colour = "Read Depth")
ggsave("CLR_depth_PCoA.png", scale = 2)
cor.test(ord_df$Axis.1, ord_df$read_depth)
# Pearson's product-moment correlation
# data:  ord_df$Axis.1 and ord_df$read_depth
# t = -6.0591, df = 28, p-value = 1.561e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.8758170 -0.5392739
# sample estimates:
#       cor 
# -0.7532072

# Perform a PCA for comparison
# Extract CLR-transformed counts as a matrix (samples x taxa)
clr_matrix <- as(otu_table(clr_phy), "matrix")
if (taxa_are_rows(clr_phy)) {
  clr_matrix <- t(clr_matrix)  # Ensure rows = samples, columns = taxa
}

# Perform PCA
pca_result <- prcomp(clr_matrix, scale. = FALSE)

# Build PCA data frame for plotting
pca_df <- as.data.frame(pca_result$x)
pca_df$SampleID <- rownames(pca_df)

# Add metadata
metadata <- as(sample_data(clr_phy), "data.frame")
metadata$SampleID <- rownames(metadata)
pca_df <- merge(pca_df, metadata, by = "SampleID")

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = read_depth)) +
  geom_point(size = 3) +
  scale_colour_viridis_c() +
  labs(
    title = "PCA of CLR-transformed Data, coloured for read depth of samples",
    x = paste0("PC1 (", round(100 * summary(pca_result)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_result)$importance[2, 2], 1), "%)")
  )
# PCoA coordinates (from your previous ordination)
pcoa_coords <- ord_df[, c("Axis.1", "Axis.2")]
rownames(pcoa_coords) <- rownames(ord_df)

# PCA coordinates
pca_coords <- pca_df[, c("PC1", "PC2")]
rownames(pca_coords) <- pca_df$SampleID

# Merge for comparison
comparison <- merge(
  pcoa_coords, pca_coords,
  by = "row.names", all = FALSE
)
colnames(comparison) <- c("SampleID", "PCoA1", "PCoA2", "PCA1", "PCA2")

# Plot comparison
ggplot(comparison, aes(x = PCoA1, y = PCA1)) +
  geom_point() +
  labs(title = "Comparison of PCoA1 and PCA1 (CLR-transformed data)")
# Yes - PCA and PCoA produce exactly the same results here.

# Final CLR plot then:
# Transform to relative abundances
rel_abund <- microbiome::transform(sslog, "compositional")
# Apply CLR transformation
clr_phy <- microbiome::transform(rel_abund, "clr")
# calculate Euclidean distances
clr_dist <- dist(t(otu_table(clr_phy)))
# Extract CLR-transformed counts and make sure rows = samples
clr_matrix <- as(otu_table(clr_phy), "matrix")
if (taxa_are_rows(clr_phy)) {
  clr_matrix <- t(clr_matrix)  # transpose so rows = samples
}
# Perform PCA
pca_result <- prcomp(clr_matrix, scale. = FALSE)

# Build PCA data frame for plotting
pca_df <- as.data.frame(pca_result$x)
pca_df$SampleID <- rownames(pca_df)

# Add metadata
metadata <- as(sample_data(clr_phy), "data.frame")
metadata$SampleID <- rownames(metadata)
pca_df <- merge(pca_df, metadata, by = "SampleID")

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Site)) +
  geom_point(aes(size = read_depth), alpha = 0.8) +
  scale_size_continuous(name = "Read depth") +
  scale_color_igv() +
  labs(title = "PCA of CLR-transformed semen data",
       subtitle = "Both Site and Read-depth are significant predictors of microbial composition.",
       x = paste0("PC1 (", round(100 * summary(pca_result)$importance[2, 1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca_result)$importance[2, 2], 1), "%)")) +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
         plot.title = element_text(face = "bold", size = 14),
         plot.margin = ggplot2::margin(1, 1, 0.5, 1, "cm"),
         plot.subtitle = element_text(size = 11),
         axis.text = element_text(size = 10, angle = 90),
         legend.text = element_text(size = 10)) +
  stat_ellipse()
ggsave("CLR_site_depth.png", width = 8, height = 6, dpi = 300)
ggsave("CLR_site_depth.pdf", width = 8, height = 6)
# Test association with Site and Line and cofounders
adonis2(clr_dist ~ Site + read_depth, data = as(sample_data(clr_phy), "data.frame"),
        by = "margin")
# Let's expand the model.
# First we should check for colinearity between variables
cor.test(pca_df$Mojo_concentration, pca_df$Mojo_Total_Motility)
# p = 9.269e-05. These variables are highly correlated, so don't use both.
cor.test(pca_df$read_depth, pca_df$mapped_boar)
# again, as expected, p=1.158e-06, highly correlated so don't use both.
cor.test(pca_df$Alive, pca_df$Mojo_Total_Motility)
# On the verge of statistical significance here - as expected.
# Now we know what to avoid in the statistical model
adonis2(clr_dist ~ Site + Line + read_depth + Mojo_concentration + Alive + Total_Motility + DFI + Acrosome_damaged,
        data = as(sample_data(clr_phy), "data.frame"),
        by = "margin")
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999

# adonis2(formula = clr_dist ~ Site + Line + read_depth + Mojo_concentration + Alive + Total_Motility + DFI + Acrosome_damaged, data = as(sample_data(clr_phy), "data.frame"), by = "margin")
# Df SumOfSqs      R2      F Pr(>F)  
# Site                1    419.4 0.04082 1.2965  0.042 *
# Line                3   1064.5 0.10360 1.0969  0.149  
# read_depth          1    465.1 0.04526 1.4377  0.012 *
# Mojo_concentration  1    338.7 0.03296 1.0471  0.328  
# Alive               1    348.8 0.03394 1.0781  0.268  
# Total_Motility      1    371.4 0.03614 1.1480  0.140  
# DFI                 1    400.3 0.03896 1.2374  0.066 .
# Acrosome_damaged    1    340.7 0.03316 1.0533  0.309  
# Residual           19   6146.1 0.59816                
# Total              29  10275.0 1.00000                
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Based on this - the most significant aspects to focus the model on are: Site and read depth.
# Hence:
adonis2(clr_dist ~ Site + read_depth,
        data = as(sample_data(clr_phy), "data.frame"),
        by = "margin")
# Here, Site explains 4.973% of the variance at p=0.006 **
# Whilst Read Depth explains 6.563% of the variance at p=0.001 ***

######################################################################
#
# Section 7: Section Bar plots and other taxonomic summaries
#
######################################################################

# We are using the dataset for all samples, except low depth (<2000 reads)
# genus_physeq_no_outliers
ntaxa(genus_physeq_no_outliers)
# 585
# How many genera would be present after filtering?
length(get_taxa_unique(genus_physeq_no_outliers, taxonomic.rank = "Genus"))
## [1] 194

# Let's try the phyla prevalence plot again on the phyloseq object with outliers removed.
# Result: a vector with one value per ASV.
prevdf3 = apply(X = otu_table(genus_physeq_no_outliers),
                MARGIN = ifelse(taxa_are_rows(genus_physeq_no_outliers), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Now you create a data.frame where each row is an ASV and columns include:
# Prevalence = number of samples where ASV is present
# TotalAbundance = total read count across all samples
# tax_table(ps0) = adds taxonomy data for each ASV
# Add taxonomy and total read counts to this data.frame
prevdf3 = data.frame(Prevalence = prevdf3,
                     TotalAbundance = taxa_sums(genus_physeq_no_outliers),
                     tax_table(genus_physeq_no_outliers))
# Compute the total and average prevalences of the features in each phyla
plyr::ddply(prevdf3, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Subset to the remaining phyla
prevdf4 = subset(prevdf3, Phylum %in% get_taxa_unique(genus_physeq_no_outliers, "Phylum"))
ggplot(prevdf4, aes(TotalAbundance, Prevalence / nsamples(genus_physeq_no_outliers),color=Phylum)) +
  # Include a guess for parameter
  # geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + 
  geom_point(size = 1, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) +
  labs(title = "Prevalence versus total counts of ASVs in Genus Amplicon dataset",
       subtitle = "Dataset has been filtered of contaminants and low prevalence ASVs") +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        axis.text.x = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "none")
ggsave("Phyla_ASV_prevalence_no_outliers.pdf", width = 8, height = 6, units = "in")

genus_glom = phyloseq::tax_glom(genus_physeq_no_outliers, "Genus", NArm = TRUE)
# Finally let's repeat the phyla prevalence plot with the agglomerated phyloseq object
prevdf5 = apply(X = otu_table(genus_glom),
                MARGIN = ifelse(taxa_are_rows(genus_glom), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Now you create a data.frame where each row is an ASV and columns include:
# Prevalence = number of samples where ASV is present
# TotalAbundance = total read count across all samples
# tax_table(ps0) = adds taxonomy data for each ASV
# Add taxonomy and total read counts to this data.frame
prevdf5 = data.frame(Prevalence = prevdf5,
                     TotalAbundance = taxa_sums(genus_glom),
                     tax_table(genus_glom))
# Compute the total and average prevalences of the features in each phyla
plyr::ddply(prevdf5, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Subset to the remaining phyla
prevdf6 = subset(prevdf5, Phylum %in% get_taxa_unique(genus_glom, "Phylum"))
ggplot(prevdf6, aes(TotalAbundance, Prevalence / nsamples(genus_glom),color=Phylum)) +
  # Include a guess for parameter
  # geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + 
  geom_point(size = 1, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) +
  labs(title = "Prevalence versus total counts of ASVs in Genus Amplicon dataset",
       subtitle = "Dataset has been filtered of contaminants and low prevalence ASVs") +
  theme_bw(base_size = 10) +
  theme(text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        axis.text.x = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "none")
# Not sure what this tells us - obviously fewer ASVS as we already knew. Only point of interest to me
# Is that the campylobacter phyla, which looked fairly populous is now formed of just two ASVs after
# genus agglomeration. Campylobacter is obviously the most prevalent and it has risen up to show it is
# present in almost half of the samples. V interesting. It is common in GI tract of animals, but has
# also been found in reproductive tract of animals and is considered pathogenic.

# Going to use the genus_physeq_no_outliers object, but transformed to relative abundance here
genus_rel=transform_sample_counts(genus_physeq_no_outliers, function(x){x / sum(x)})
# agglomerate taxa
phyglom <- tax_glom(genus_rel, taxrank = 'Phylum', NArm = FALSE)
phymelt <- psmelt(phyglom)
# change to character for easy-adjusted level
phymelt$phylum <- as.character(phymelt$Phylum)

# Let's just get phyla stats here for semen:
# Filter to just Firmicutes and desired sample type
firm_df <- phymelt |>
  filter(Phylum == "Firmicutes", sampletype == "semen")

# Calculate mean and SD of Firmicutes relative abundance within semen samples
firm_df |>
  group_by(sampletype) |>
  summarise(mean_abundance = mean(Abundance),
            sd_abundance = sd(Abundance),
            n = n())
# semen, mean_abundance firmicutes 0.772, sd 0.176

# First step is to try a faceted plot of all sample types
stphymelt <- phymelt |>
  group_by(Sample, Phylum) |>
  mutate(median=median(Abundance))
#to get the same rows together
stphymelt_sum <- stphymelt |>
  group_by(Sample, sampletype, Phylum) |>
  summarise(Abundance=sum(Abundance))
# Example: order by Firmicutes abundance
ordered_samples <- stphymelt_sum |>
  filter(Phylum == "Firmicutes") |>
  arrange(sampletype, desc(Abundance)) |>
  pull(Sample)
stphymelt_sum$Sample <- factor(stphymelt_sum$Sample, levels = ordered_samples)
# Now make the plot
stphymelt_sum |>
  group_by(Phylum) |>
  mutate(mean_abund = mean(Abundance)) |>
  ungroup() |>
  mutate(Phylum = fct_reorder(Phylum, mean_abund)) |>
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1, aes(fill=Phylum)) + 
  labs(x = "", y = "Relative Abundance (%)") +
  facet_grid(~ sampletype, scales= "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_igv() +
    theme_bw(base_size = 10) +
    theme(strip.background = element_blank(),
          strip.text.x.top = element_text(size = 10, angle = 90),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          axis.text.x.bottom = element_blank(),
          text = element_text(size = 10),
          plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10)) +
    labs(title = "Relative Abundance of Phyla in Boar Dataset",
       subtitle = "Samples have been grouped by Sample Type")
ggsave("Phyla_sampletype.png", width = 7, height = 5, dpi = 300)
ggsave("Phyla_sampletype.pdf", width = 8, height = 6, units = "in")
# Now let's look at the semen samples by Site
phymelt <- psmelt(phyglom)
# change to character for easy-adjusted level
phymelt$phylum <- as.character(phymelt$Phylum)
phymelt <- phymelt |>
  filter(sampletype == "semen") |>
  group_by(Sample, Phylum) |>
  mutate(median=median(Abundance))
# select group median > 1
keep <- unique(phymelt$Phylum[phymelt$median > 0.01])
phymelt$Phylum[!(phymelt$Phylum %in% keep)] <- "< 1%"
#to get the same rows together
phymelt_sum <- phymelt |>
  filter(sampletype == "semen") |>
  group_by(Sample, Site, Phylum) |>
  summarise(Abundance=sum(Abundance))
# Example: order by Firmicutes abundance
ordered_samples <- phymelt_sum |>
  filter(Phylum == "Firmicutes") |>
  arrange(Site, desc(Abundance)) |>
  pull(Sample)
phymelt_sum$Sample <- factor(phymelt_sum$Sample, levels = ordered_samples)
phymelt_sum |>
  group_by(Phylum) |>
  mutate(mean_abund = mean(Abundance)) |>
  ungroup() |>
  mutate(Phylum = fct_reorder(Phylum, mean_abund)) |>
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1, aes(fill=Phylum)) + 
  labs(x="", y="%") +
  facet_wrap(~ Site, scales= "free_x", nrow=1) +
  theme_classic(base_size = 10) + 
  scale_fill_igv() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(strip.background = element_blank(),
        strip.text.x.top = element_text(size = 10),
        strip.text = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text.x.bottom = element_text(angle = -90)) +
  labs(title = "Relative Abundance of Phyla in Semen Samples",
       subtitle = "Outlier samples are restricted to Boston",
       caption = "Sample 0839K is outlier in beta-diversity plots",
       y = "Relative Abundance (%)")
ggsave("Grouped_semen_phyla_bar.pdf", width = 8, height = 6, units = "in")

# Genus for semen
genusglom <- tax_glom(genus_rel, taxrank = 'Genus', NArm = FALSE)
ps.melt.genus <- psmelt(genusglom)
# change to character for easy-adjusted level
ps.melt.genus$Genus <- as.character(ps.melt.genus$Genus)
ps.melt.genus <- ps.melt.genus |>
  filter(sampletype == "semen") |>
  group_by(Sample, Site, Class, Order, Family, Genus) |>
  mutate(median=median(Abundance))
# Create a new label column
ps.melt.genus$Genus_label <- with(ps.melt.genus,
                                      ifelse(
                                        is.na(Genus) | Genus == "uncultured" | Genus == "" | Genus == "NA",
                                        ifelse(!is.na(Family) & Family != "",
                                               paste0("Family: ", Family),
                                               ifelse(!is.na(Order) & Order != "",
                                                      paste0("Order: ", Order),
                                                      ifelse(!is.na(Class) & Class != "",
                                                             paste0("Class: ", Class),
                                                             "Unclassified"))),
                                        Genus))
# select group median > 1
keep <- unique(ps.melt.genus$Genus_label[ps.melt.genus$median > 0.075])
ps.melt.genus$Genus_label[!(ps.melt.genus$Genus_label %in% keep)] <- "< 7.5%"
#to get the same rows together
ps.melt.genus_sum <- ps.melt.genus |>
  group_by(Sample,Site,Genus_label) |>
  summarise(Abundance=sum(Abundance))
# Example: order by Firmicutes abundance
ordered_samples <- ps.melt.genus_sum |>
  filter(Genus_label == "< 7.5%") |>
  arrange(Site, desc(Abundance)) |>
  pull(Sample)
ps.melt.genus_sum$Sample <- factor(ps.melt.genus_sum$Sample, levels = ordered_samples)
ps.melt.genus_sum |>
  group_by(Genus_label) |>
  mutate(mean_abund = mean(Abundance)) |>
  ungroup() |>
  mutate(Genus_label = fct_reorder(Genus_label, mean_abund)) |>
  ggplot(aes(x = Sample, y = Abundance, fill = Genus_label)) + 
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1, aes(fill=Genus_label),
           show.legend = F) + 
  facet_wrap(~ Site, scales= "free_x", nrow=1) +
  theme_classic() + 
  scale_fill_igv() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_blank(),
        strip.text.x.top = element_text(size = 8),
        strip.text = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10)) +
  labs(title = "Relative Abundance of Semen Genera",
       subtitle = "Rare taxa with median relative abundance < 7.5% have been agglomerated",
       x = "Samples",
       y = "Relative Abundance (%)")
ggsave("semen_genus_bar.pdf", width = 8, height = 6, units = "in")
# Issue here is that the samples are too full of variation - and difficult to visualise main taxa.

# Reset the objects
genusglom <- tax_glom(genus_rel, taxrank = 'Genus', NArm = FALSE)
ps.melt.genus <- psmelt(genusglom)
#to get the same rows together
ps.melt.genus_sum <- ps.melt.genus |>
  group_by(Sample,Site, Class, Order, Family, Genus) |>
  summarise(Abundance=sum(Abundance))
# Create a new label column
ps.melt.genus_sum$Genus_label <- with(ps.melt.genus_sum,
                                      ifelse(
                                        is.na(Genus) | Genus == "uncultured" | Genus == "" | Genus == "NA",
                                        ifelse(!is.na(Family) & Family != "",
                                               paste0("Family: ", Family),
                                               ifelse(!is.na(Order) & Order != "",
                                                      paste0("Order: ", Order),
                                                      ifelse(!is.na(Class) & Class != "",
                                                             paste0("Class: ", Class),
                                                             "Unclassified"))),
                                        Genus))
# Calculate mean abundance per informative label
top_genera <- ps.melt.genus_sum %>%
  group_by(Genus_label) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 10) %>%
  pull(Genus_label)
# Set all other genera to "Other"
ps.melt.genus_sum$Genus_label[!(ps.melt.genus_sum$Genus_label %in% top_genera)] <- "Other"
# After assigning "Other" to low-abundance genera:
ps.melt.genus_sum$Genus_label <- as.character(ps.melt.genus_sum$Genus_label)
# Set factor levels for plotting
ps.melt.genus_sum$Genus_label <- factor(
  ps.melt.genus_sum$Genus_label,
  levels = c(top_genera, "Other")
)
# Plot
ggplot(ps.melt.genus_sum, aes(x = Sample, y = Abundance, fill = Genus_label)) +
  geom_bar(stat = "identity", colour = "black", width = 0.7) +
  facet_wrap(~ Site, scales= "free_x", nrow=1) +
  theme_classic() +
  scale_fill_igv() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  labs(title = "Top 10 Genera in Semen Samples",
       subtitle = "All other genera grouped as 'Other'",
       caption = "Rectangles within 'Other' represent distinct taxa")
ggsave("Semen_top_10.png", scale = 2)
# Notes on the top 10 genera
# E coli - the most abundant taxa in the entire dataset (not just semen), this is likely because
# 1: it was present in mock samples, and very overrepresented in cell mocks. Plus it was also
# present in negative control - so some wariness about this ASV - difficult to interpret.
# We know that this is biologically plausible as it is easy to culture from semen samples yet is
# also a potential contaminant from gut/faecal and there is abundant presence in other niches .

# Family: Carnobacteriaceae, Genus unclassified. Most abundant in semen samples. 2nd in all data.
# BLAST search reveals uncultured organisms as top hits.
# Carnobacteriaceae, a family of the order Lactobacillales, within the phylum Firmicutes, embraces 
# the genera Carnobacterium, Alkalibacterium, Allofustis, Alloiococcus, Atopobacter, Atopococcus, 
# Atopostipes, Bavariicoccus, Desemzia, Dolosigranulum, Granulicatella, Isobaculum, Jeotgalibaca, 
# Lacticigenium, Marinilactibacillus, Pisciglobus, and Trichococcus.
# Gram positive and facultative anaerobes.
# These have been mainly characterised via sequencing - so this result is not surprising. Prior to
# NGS approaches - remarked upon as abberant lactobacilli. Common in animal hosts, and as other
# lactic acid bacteria - main substrate is glucose. Definitely found in other boar semen microbiome
# studies - eg Ngo et al 2023.
# There is not much to go on here.

# Dolosicoccus - Gram postive facultative anaerobe. relatively newly characterised - offshoot of 
# Streptococci (i.e. as taxonomy has expanded alongside NGS approaches). Closely related to other 
# members of this family Aerococcaceae. Dolosicoccus has not been reported in previous boar semen 
# studies, although there have been examples of uncharacterised members of this family eg: Ngo et al 2023. 
# In addition it has been reported in boar intestinal microbiota at high prevalence McCormack et al 2019.

# Streptococcus - commonly found as part of animal and human seminal microbiota. Despite this found
# in other niches here and potential for contaminant.

# Class: Bacilli - unclassified, order, family and genera. BLAST search reveals many uncultured rRNA 
# sequences and a few Aerococcus sprinkled into the results. Appears, as with Dolosicoccus, related
# but not enough sequence difference to annotate at higher levels. However, perfectly reasonable
# that more bacilli are in the microbiota, as these are well reported in previous studies.

# W5303 - not much information about this mysterious bacteria out there, except that it regularly
# appears in sequencing data from animal reproductive microbiomes - particularly semen, including,
# cattle (Cojkic 2021), horses (Malaluang 2024) and ewes. Family XI is a placeholder name, representing
# a clade of gram positive anaerobes - previously listed in the order of Clostridiales. Other closely
# related bacteria include Finegoldia and Anaerococcus which are all well represented in NGS
# output of animal related reproductive microbiomes. These tend to be gram positive obligate anaerobes.

# Lactobacillis - well characterised and expected here.

# Clostridium sensus stricto 1. This is potentially an environmental contaminant, as it appeared 
# in high prevalence on gloves and skin swabs. As clostridia produce spores, these will be environmentally
# resistant and contamination will be easy.Action here - review with venn diagram, and possibly look
# at removing this from semen data.

# Psychrobacter - an obvious contender for contaminant organism. Environmentally resistant and present
# in skin and glove swabs.

# Gram positive and aerobic - commonly found in semen microbiome studies. However, still found in 
# skin and glove swabs. Despite plausibility could very well be introduced environmentally.

# Let's try a bubble plot
library(forcats)
# Reset the objects
genusglom <- tax_glom(genus_rel, taxrank = 'Genus', NArm = FALSE)
ps.melt.genus <- psmelt(genusglom)
#to get the same rows together
ps.melt.genus <- ps.melt.genus |>
  group_by(Sample, sampletype, Site, Class, Order, Family, Genus)
# Create a new label column
ps.melt.genus$Genus_label <- with(ps.melt.genus,
                                      ifelse(
                                        is.na(Genus) | Genus == "uncultured" | Genus == "" | Genus == "NA",
                                        ifelse(!is.na(Family) & Family != "",
                                               paste0("Family: ", Family),
                                               ifelse(!is.na(Order) & Order != "",
                                                      paste0("Order: ", Order),
                                                      ifelse(!is.na(Class) & Class != "",
                                                             paste0("Class: ", Class),
                                                             "Unclassified"))),
                                        Genus))

# Prep genus table (already agglomerated and relative abundance)
top_genera <- ps.melt.genus |>
  group_by(Genus_label) |>
  summarise(mean_abund = mean(Abundance)) |>
  arrange(desc(mean_abund)) |>
  slice_head(n = 20)

bubble_data <- ps.melt.genus |>
  filter(Genus_label %in% top_genera$Genus_label)
# Count true zero abundances
sum(bubble_data$Abundance == 0) # 224
# Or better: check the minimum non-zero value
min(bubble_data$Abundance[bubble_data$Abundance > 0]) #0.00014856
# Ok, I think zeros are being plotted as bubbles here.
# Try this
bubble_data <- bubble_data |>
  mutate(present = ifelse(Abundance > 0.001, "Present", "Absent"))

# 1.1  work out the ordering vector
semen_order <- ps.melt.genus |>
  filter(sampletype == "semen") |>          # look only at semen samples
  group_by(Genus_label) |>
  summarise(semen_mean = mean(Abundance), .groups = "drop") |>
  arrange(desc(semen_mean)) |>
  pull(Genus_label)

# 1.2  convert Genus_label into an ordered factor
bubble_data <- bubble_data |>
  mutate(Genus_label = factor(Genus_label, levels = semen_order))
# Before plotting consider using raw counts
libsize <- sample_sums(genus_physeq_no_outliers)              # your counts object
depth_df <- data.frame(Sample = names(libsize),
                       TotalReads = as.numeric(libsize))
bubble_data <- bubble_data |>
  left_join(depth_df, by = "Sample") |>
  mutate(ReadCount = round(Abundance * TotalReads))

ggplot(bubble_data, aes(x = Sample, y = Genus_label)) +
  geom_point(aes(size = Abundance, fill = Abundance, shape = present), alpha = 0.8) +
  scale_shape_manual(values = c("Absent" = NA, "Present" = 21)) +
  scale_size_area(max_size = 6) +
  facet_grid(~ sampletype,
            scales = "free_x",
            space  = "free_x",
            drop   = TRUE) +
  scale_fill_viridis_c(labels = scales::percent_format(accuracy = 0.1)) +
  theme_bw(base_size = 10) +
  labs(x = "Sample",
       y = "Genus",
       title = "Bubble plot of Top 20 Genera in Boar Dataset",
       subtitle = "Bubble size and colour represents relative abundance",
       fill = "Relative Abundance",
       size = "Relative Abundance",
       shape = "Presence") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x.top = element_text(size = 8, angle = 90),
        strip.text = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text.x.bottom = element_blank())
ggsave("Top20_bubble.pdf", width = 10, height = 8, units = "in")  
ggsave("Top20_bubble.png", width = 9, height = 7, dpi = 300)
# Now let's explore dot plots
#to get the same rows together and summarise mean abundance
genus_summary <- ps.melt.genus |>
  filter(sampletype == "semen") |>
  group_by(Site, Genus_label) |>
  summarise(mean_abund = mean(Abundance), sd_abund = sd(Abundance)) |>
  ungroup() |>
  group_by(Genus_label) |>
  filter(any(mean_abund > 0.009))
# Plot
genus_summary |>
  filter(mean_abund > 0) |>  # remove zero values
  ggplot(aes(x = mean_abund, y = fct_reorder(Genus_label, mean_abund), color = Site)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = mean_abund - sd_abund, xmax = mean_abund + sd_abund),
                 position = position_dodge(width = 0.5), height = 0.2) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Mean Relative Abundance", y = "Genus",
       title = "Mean Relative Abundance of Semen Genera, by Site",
       subtitle = "Filtering threshold > 0.009 mean abundance") +
  theme_classic(base_size = 10) +
  scale_color_igv() +
  theme(legend.title = element_text(size = 10, face = "bold"),
        text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10))
ggsave("Semen_dot.pdf", width = 8, height = 6, units = "in")
# This is interesting, most abundance genera are shared - mean abundances do not differ enormously.
# Real differences occur by site at low abundance. Need to tie this in somehow with unique genera.
# try this again, but this time include all sample types
# Now let's explore dot plots
#to get the same rows together and summarise mean abundance
genus_summary <- ps.melt.genus |>
  group_by(sampletype, Genus_label) |>
  summarise(mean_abund = mean(Abundance), sd_abund = sd(Abundance)) |>
  filter(mean_abund > 0.01)  # Optional filter
# Plot
genus_summary |>
  ggplot(aes(x = mean_abund, y = fct_reorder(Genus_label, mean_abund), color = sampletype)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = mean_abund - sd_abund, xmax = mean_abund + sd_abund),
                 position = position_dodge(width = 0.5), height = 0.2) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Mean Relative Abundance", y = "Genus",
       title = "Mean Relative Abundance of Semen Genera, by Site",
       subtitle = "Filtering threshold > 0.1 mean abundance") +
  theme_classic(base_size = 10) +
  scale_color_igv() +
  theme(legend.title = element_text(size = 10, face = "bold"),
        text = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10))
ggsave("Semen_dot.pdf", width = 8, height = 6, units = "in")

# Finally, genus level for the other samples
genusglom <- tax_glom(genus_rel, taxrank = 'Genus', NArm = FALSE)
ps.melt.genus <- psmelt(genusglom)
# change to character for easy-adjusted level
ps.melt.genus$Genus <- as.character(ps.melt.genus$Genus)

ps.melt.genus <- ps.melt.genus |>
  filter(sampletype != "semen") |>
  group_by(Sample, Genus) |>
  mutate(median=median(Abundance))
# select group median > 1
keep <- unique(ps.melt.genus$Genus[ps.melt.genus$median > 0.05])
ps.melt.genus$Genus[!(ps.melt.genus$Genus %in% keep)] <- "< 5%"
#to get the same rows together
ps.melt.genus_sum <- ps.melt.genus |>
  group_by(Sample,sampletype,Genus) |>
  summarise(Abundance=sum(Abundance))
ps.melt.genus_sum |>
  group_by(Genus) |>
  mutate(mean_abund = mean(Abundance)) |>
  ungroup() |>
  mutate(Genus = fct_reorder(Genus, mean_abund)) |>
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", colour = "black", width = 0.7, aes(fill=Genus)) + 
  labs(x="", y="%") +
  facet_wrap(~ sampletype, scales= "free_x", nrow=1) +
  theme_classic() + 
  scale_fill_igv() +
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90)) +
  labs(title = "Relative Abundance of Genera in Associated Samples",
       subtitle = "Rare taxa with median relative abundance < 5% have been agglomerated")
ggsave("other_genus_bar.png", scale = 2)

##################################################################################
#
# 8. Exploration of shared taxa across niches
#
##################################################################################
# Begin this section by exploring using a simple Venn plot.

# Agglomerate at genus level (or other, e.g., "Phylum", "ASV" for OTU level)
phy_genus <- tax_glom(genus_rel, taxrank = "Genus")
# Extract presence/absence per sample
otu_mat <- otu_table(phy_genus)
otu_mat[otu_mat > 0] <- 1  # Convert to presence/absence
otu_mat <- as.data.frame(t(otu_mat))  # Samples as rows
sample_data_df <- as(sample_data(phy_genus), "data.frame")

# Add sample type info
otu_mat$sampletype <- sample_data_df$sampletype
# Split by sampletype and get union of present taxa per group
venn_input <- otu_mat %>%
  group_by(sampletype) %>%
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames("sampletype") %>%
  as.matrix()

# Convert sums to presence/absence (≥1 means present in that type)
venn_input[venn_input >= 1] <- 1

# Create list for Venn diagram
venn_list <- apply(venn_input, 1, function(x) names(which(x == 1)))

## Current names (check with):
names(venn_list)
# [1] "boarswab"  "extender"  "gloveswab" "semen" 

# --- give them whatever labels you prefer -----------------------------
names(venn_list) <- c("Boar skin swab",
                      "Extender",
                      "Glove swab",
                      "Semen")


# Plot Venn
ggvenn(venn_list, text_size = 3, set_name_size = 3, stroke_size   = 0.7) +
  labs(title = "Venn diagram of shared taxa within 4 niches at 2 boar studs",
       subtitle = "A majority of taxa are unique to semen") +
  theme_classic(base_size = 10) +
  scale_color_igv() +
  theme(text = element_text(size = 10),
        plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10))
ggsave("niche_venn.png", width = 7, height = 5, dpi = 300)
ggsave("niche_venn.pdf", width = 8, height = 6, units = "in")
# The Venn diagram shows us that we clearly have shared taxa and some unique to semen,
# but we need to explore more about these particular taxa.

# First look at abundance bar plots of shared taxa.

# Find genera shared by all groups (intersection)
shared_genera <- Reduce(intersect, venn_list)
print(shared_genera)
# Find genera unique to a group (e.g., "Semen")
unique_semen <- setdiff(venn_list[["Semen"]], 
                        unlist(venn_list[names(venn_list) != "Semen"]))
print(unique_semen)

# 1. Agglomerate at genus level and melt
ps.melt.genus <- psmelt(genusglom)
ps.melt.genus$Genus <- as.character(ps.melt.genus$Genus)
ps.melt.genus$Best_Taxonomy <- ifelse(
  !is.na(ps.melt.genus$Genus) & ps.melt.genus$Genus != "",
  paste0(ps.melt.genus$Genus),
  ifelse(
    !is.na(ps.melt.genus$Family) & ps.melt.genus$Family != "",
    paste0(ps.melt.genus$Family),
    ifelse(
      !is.na(ps.melt.genus$Order) & ps.melt.genus$Order != "",
      paste0("Order: ", ps.melt.genus$Order),
      ifelse(
        !is.na(ps.melt.genus$Class) & ps.melt.genus$Class != "",
        paste0("Class: ", ps.melt.genus$Class),
        "Unclassified"
        )
      )
    )
  )
ps.melt.genus$Best_Taxonomy <- as.character(ps.melt.genus$Best_Taxonomy)
# 2. Calculate mean abundance per Genus per sampletype
genus_abund <- ps.melt.genus %>%
  group_by(sampletype, Best_Taxonomy) %>%
  summarise(mean_abund = mean(Abundance),sd_abund = sd(Abundance), .groups = "drop")

# 3. Identify rare genera (mean abundance <1% across all sampletypes)
rare_genera <- genus_abund |>
  group_by(Best_Taxonomy) |>
  summarise(overall_mean = mean(mean_abund)) |>
  filter(overall_mean < 0.01) |>
  pull(Best_Taxonomy)

# 4. Assign "< 1%" to rare genera
genus_abund$Best_Taxonomy[genus_abund$Best_Taxonomy %in% rare_genera] <- "< 1%"

# 5. Remove "< 1%" from further analysis (optional, for clarity)
genus_abund_filtered <- genus_abund |>
  filter(Best_Taxonomy != "< 1%")

# 6. Presence/absence table for shared genera logic
genus_pa <- ps.melt.genus %>%
  group_by(sampletype, Best_Taxonomy) %>%
  summarise(present = sum(Abundance) > 0, .groups = "drop") %>%
  filter(present)
# 7. Find genera present in semen and at least one other sampletype
genus_in_semen <- genus_pa |>
  group_by(Best_Taxonomy) |>
  summarise(
    in_semen = any(sampletype == "semen"),
    n_sampletypes = n_distinct(sampletype)
  ) |>
  filter(in_semen, n_sampletypes >=2) |>
  pull(Best_Taxonomy)

# 8. Now filter the abundance table for genera of interest
shared_with_semen_abund <- genus_abund_filtered |>
  filter(Best_Taxonomy %in% genus_in_semen)
# 9. Order Genus factor by overall mean abundance
genus_order <- shared_with_semen_abund %>%
  group_by(Best_Taxonomy) %>%
  summarise(total = sum(mean_abund)) %>%
  arrange(desc(total)) %>%
  pull(Best_Taxonomy)
shared_with_semen_abund$Best_Taxonomy <- factor(shared_with_semen_abund$Best_Taxonomy, levels = genus_order)
# 10. Plot
ggplot(shared_with_semen_abund, aes(x = sampletype, y = mean_abund, fill = Best_Taxonomy)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Best_Taxonomy, scales = "free_y") +
  labs(title = "Abundance of Genera shared between Semen and at least one other sample type",
       subtitle = "Genera with mean abundance per sample type <1% have been filtered out",
       x = "Sample Type",
       y = "Mean relative abundance",
       fill = "Genus") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_igv()
ggsave("Abundance_shared_genera.png", scale = 2)
# Issue with this plot is that information on directionality is unavailable.
# We might guess that taxa that are more abundant in semen are more likely to originate
# in these samples, but this is still just an assumption.
# I have made notes on most of these, but in addition we can see that
# Jeotgalicoccus - not associated with host microbiomes - potentially
# environmental
# Acidothermus - potentially environmental
# HT002 - a lactobacillus species
# Turicibacter - common gut commensal of animals - found in faeces, so 
# lots of potential for environmental contamination.
# Terrisporobacter - gut commensal and produces endospores - likely
# environmental contaminant.
# Jeotgalibaca - isolated from pig joints! Full habitat unknown.
# Enhydrobacter - likely contaminant as predominantly in extender.
# Now we explore taxa Unique to Semen
# 1. Calculate mean abundance per Genus per sampletype
genus_abund <- ps.melt.genus |>
  group_by(sampletype, Best_Taxonomy) |>
  summarise(mean_abund = mean(Abundance), sd_abund = sd(Abundance), .groups = "drop")
# 2. Presence/absence table
genus_pa <- ps.melt.genus |>
  group_by(sampletype, Best_Taxonomy) |>
  summarise(present = sum(Abundance) > 0, .groups = "drop") |>
  filter(present)
# 3. Find genera unique to semen
genus_unique_semen <- genus_pa |>
  group_by(Best_Taxonomy) |>
  summarise(
    unique_semen = all(sampletype == "semen"),
    n_sampletypes = n_distinct(sampletype)
  ) |>
  filter(unique_semen, n_sampletypes == 1) |>
  pull(Best_Taxonomy)
# 4. Calculate prevalence in semen samples
prevalence <- ps.melt.genus |>
  filter(sampletype == "semen") |>
  group_by(Best_Taxonomy) |>
  summarise(prevalence = sum(Abundance > 0), .groups = "drop")
# 5. Filter abundance table for unique semen genera, only for semen sampletype
semen_unique_abund <- genus_abund |>
  filter(Best_Taxonomy %in% genus_unique_semen, sampletype == "semen") |>
  left_join(prevalence, by = "Best_Taxonomy")
# 6. Filter for genera present in >3 semen samples
semen_unique_abund_filtered <- semen_unique_abund |>
  filter(prevalence > 2)
# 7. Plots
semen_unique_abund_filtered |>
  ggplot(aes(y = prevalence, x = mean_abund, label = Best_Taxonomy)) +
  geom_point() +
  geom_text_repel(max.overlaps = 20) +
  labs(title = "Unique-to-Semen Genera: Prevalence vs Mean Abundance",
       subtitle = "Mean Abundance calculated across all semen samples",
       x = "Prevalence (# of semen samples)",
       y = "Mean Abundance") +
  theme_bw()
ggsave("semen_unique_genera_repel.png", scale = 3)
# Also write the semen_unique_abundance table to csv
write.csv(semen_unique_abund, file = "semen_unique_abundance.csv")

# Recalculate abundance, but this time per sample, in order to capture 
# rare taxa for comparison
# 1. Calculate mean abundance per Genus per Sample for those unique to semen,
semen_unique_sample_abund <- ps.melt.genus |>
  filter(sampletype == "semen", Best_Taxonomy %in% genus_unique_semen) |>
  group_by(Sample, Best_Taxonomy) |>
  summarise(abund = sum(Abundance), .groups = "drop")
# 2. Wide format for heatmap
abund_matrix <- semen_unique_sample_abund |>
  tidyr::pivot_wider(names_from = Best_Taxonomy, values_from = abund, values_fill = 0) |>
  column_to_rownames("Sample")
# Ensure rownames of annotation match rownames of your abundance matrix
annotation_row <- Metadata_semen |>
  filter(sampleid %in% rownames(abund_matrix)) |>
  column_to_rownames("sampleid") |>
  select(Site, Line, depth_binned)  # or whatever columns you want to annotate

log_abund_matrix <- log10(t(as.matrix(abund_matrix)) + 1e-5)  # Add small value to avoid log(0)
# Fix legend scaling
breaks <- seq(-5, -1, by = 1)
labels <- c("1e-5", "1e-4", "1e-3", "1e-2", "1e-1")
annotation_row$depth_binned <- factor(annotation_row$depth_binned,
                                      levels = c("(2e+03,5e+03]", "(5e+03,1e+04]", "(1e+04,6e+04]"),
                                      labels = c("Low", "Med", "High"))
colnames(annotation_row)[colnames(annotation_row) == "depth_binned"] <- "Read Depth"
pheat <- pheatmap(
  log_abund_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_row,
  colour = viridis(100),
  fontsize_row = 6,
  angle_col = 45,
  cutree_cols = 6,
  annotation_colors = list(
    Site = c("Boston" = "steelblue", "Willingham" = "firebrick"),
    Line = c("L03_LW" = "gold", "L04_LR" = "purple", "L19_WD" = "green", "327_Hamp" = "orange")
  ),
  legend_breaks = breaks,
  legend_labels = labels,
  main = "Log10 Abundance of Unique-to-Semen Genera by Sample"
)
ggsave("heatmap_unique_semen.png", plot = pheat, width = 10, height = 15, dpi = 300)
ggsave("heatmap_unique_semen.pdf", plot = pheat, width = 10, height = 15)
# Issue with this plot is that genera are all extremely low abundance.
# hence log10 transformation.
# These data may highlight v low biomass signal. Or it could reveal
# sequencing artefacts. Or both - which, lets face it, is extremely likely.

# Let's try the same plot but at higher taxonomic levels, to see if there are any patterns or
# particular groups we can see.

# Use the ps.melt.genus object, but create a new best taxonomy at higher level label to plot with.
# Create a best taxonomy label to include and annotate those ASVs annotated at 
# less informative levels
ps.melt.genus$Best_Family_Taxonomy <- ifelse(
    !is.na(ps.melt.genus$Family) & ps.melt.genus$Family != "",
    paste0(ps.melt.genus$Family),
    ifelse(
      !is.na(ps.melt.genus$Order) & ps.melt.genus$Order != "",
      paste0("Order: ", ps.melt.genus$Order),
      ifelse(
        !is.na(ps.melt.genus$Class) & ps.melt.genus$Class != "",
        paste0("Class: ", ps.melt.genus$Class),
        "Unclassified"
      )
    )
  )
ps.melt.genus$Best_Family_Taxonomy <- as.character(ps.melt.genus$Best_Family_Taxonomy)
# 1. Calculate mean abundance per Genus per Sample for those unique to semen,
semen_unique_sample_abund_family <- ps.melt.genus |>
  filter(sampletype == "semen", Best_Taxonomy %in% genus_unique_semen) |>
  group_by(Sample, Best_Family_Taxonomy) |>
  summarise(abund = sum(Abundance), .groups = "drop")
# 2. Wide format for heatmap
abund_matrix_family <- semen_unique_sample_abund_family |>
  tidyr::pivot_wider(names_from = Best_Family_Taxonomy, values_from = abund, values_fill = 0) |>
  column_to_rownames("Sample")

log_abund_matrix_family <- log10(t(as.matrix(abund_matrix_family)) + 1e-5)  # Add small value to avoid log(0)

fampheat <- pheatmap(
  log_abund_matrix_family,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_row,
  colour = viridis(100),
  fontsize_row = 7,
  fontsize_col = 5,
  angle_col = 45,
  cutree_cols = 4,
  annotation_colors = list(
    Site = c("Boston" = "steelblue", "Willingham" = "firebrick"),
    Line = c("L03_LW" = "gold", "L04_LR" = "purple", "L19_WD" = "green", "327_Hamp" = "orange")
  ),
  legend_breaks = breaks,
  legend_labels = labels,
  main = "Log10 Abundance of Unique-to-Semen Bacterial Families"
)
ggsave("family_unique_heat.pdf", plot = fampheat, width = 8, height = 10, units = "in")
# Let's try order, just in case.
# Create a best taxonomy label to include and annotate those ASVs annotated at 
# less informative levels
ps.melt.genus$Best_Order_Taxonomy <- ifelse(
    !is.na(ps.melt.genus$Order) & ps.melt.genus$Order != "",
    paste0(ps.melt.genus$Order),
    ifelse(
      !is.na(ps.melt.genus$Class) & ps.melt.genus$Class != "",
      paste0("Class: ", ps.melt.genus$Class),
      "Unclassified"
    )
  )
ps.melt.genus$Best_Order_Taxonomy <- as.character(ps.melt.genus$Best_Order_Taxonomy)
# 1. Calculate mean abundance per Genus per Sample for those unique to semen,
semen_unique_sample_abund_order <- ps.melt.genus |>
  filter(sampletype == "semen", Best_Taxonomy %in% genus_unique_semen) |>
  group_by(Sample, Best_Order_Taxonomy) |>
  summarise(abund = sum(Abundance), .groups = "drop")
# 2. Wide format for heatmap
abund_matrix_order <- semen_unique_sample_abund_order |>
  tidyr::pivot_wider(names_from = Best_Order_Taxonomy, values_from = abund, values_fill = 0) |>
  column_to_rownames("Sample")

log_abund_matrix_order <- log10(t(as.matrix(abund_matrix_order)) + 1e-5)  # Add small value to avoid log(0)

pheat <- pheatmap(
  log_abund_matrix_order,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_row,
  colour = viridis(100),
  fontsize_row = 7,
  fontsize_col = 6,
  angle_col = 45,
  cutree_cols = 4,
  annotation_colors = list(
    Site = c("Boston" = "steelblue", "Willingham" = "firebrick"),
    Line = c("L03_LW" = "gold", "L04_LR" = "purple", "L19_WD" = "green", "327_Hamp" = "orange")
  ),
  legend_breaks = breaks,
  legend_labels = labels,
  main = "Log10 Abundance of Unique-to-Semen Microbial Orders"
)
# Again, very similar pattern observed.
ggsave("Orders_unique_heatmap.pdf", plot = pheat, width = 8, height = 10, units = "in")
# Let's look at the 4 most prevalent genera.
# Step 1: Filter to 4 target genera
target_genera <- c("Fastidiosipila", "Actinobaculum", "Macrococcus", "Providencia")
subset_abund <- ps.melt.genus |>
  filter(sampletype == "semen", Genus %in% target_genera) |>
  group_by(Sample, Genus) |>
  summarise(abund = sum(Abundance), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Genus, values_from = abund, values_fill = 0) |>
  column_to_rownames("Sample")
# Step 2: Log transform (optional, to dampen scale differences)
log_abund <- log10(subset_abund + 1e-5)
# Step 3: Calculate distance and clustering
sample_dist <- dist(log_abund, method = "euclidean")
sample_clust <- hclust(sample_dist, method = "ward.D2")
# Step 4: Dendrogram with metadata annotations (if wanted)
# E.g. ggdendro or ggtree, or use heatmap with clustering
meta_subset <- Metadata_semen |>
  filter(sampleid %in% rownames(log_abund)) |>
  column_to_rownames("sampleid") |>
  select(Site, Line, depth_binned)
meta_subset$depth_binned <- factor(meta_subset$depth_binned,
                                   levels = c("(2e+03,5e+03]", "(5e+03,1e+04]", "(1e+04,6e+04]"),
                                   labels = c("Low", "Med", "High"))
colnames(meta_subset)[colnames(meta_subset) == "depth_binned"] <- "Read Depth"
# Fix legend scaling
breaks <- seq(-5, -1, by = 1)
labels <- c("1e-5", "1e-4", "1e-3", "1e-2", "1e-1")

abund_heat <- pheatmap(
  t(log_abund),  # transpose to get taxa as rows
  cluster_cols = TRUE,
  cluster_rows = FALSE,
  annotation_col = meta_subset,
  annotation_colors = list(
    Site = c("Boston" = "steelblue", "Willingham" = "firebrick"),
    Line = c("L03_LW" = "gold", "L04_LR" = "purple", "L19_WD" = "green", "327_Hamp" = "orange")
  ),
  colour = viridis::viridis(100),
  legend_breaks = breaks,
  legend_labels = labels,
  main = "Log10 Abundance of Key Genera (Unique-to-Semen)"
)
ggsave("key_abund_heat.pdf", plot = abund_heat, width = 8, height = 6, units = "in")
# Actinobaculum shows a distinct pattern with Site and Line. Let's plot individually.
# Extract Bacilli abundance
Actinobaculum_abund <- as.data.frame(log_abund_matrix["Actinobaculum", ])
Actinobaculum_abund$Sample <- rownames(Actinobaculum_abund)
colnames(Actinobaculum_abund)[1] <- "log10_abund"

# Add metadata
Actinobaculum_abund <- left_join(Actinobaculum_abund, Metadata_semen, by = c("Sample" = "sampleid"))

# Quick comparison by Site
ggplot(Actinobaculum_abund, aes(x = Site, y = log10_abund, fill = Site)) +
  geom_boxplot(width = 0.3, show.legend = F) +
  geom_jitter(width = 0.15, show.legend = F) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_y_continuous(labels = scales::math_format(10^.x)) +
  theme_classic(base_size = 10) +
  scale_fill_igv() +
  theme(text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")) +
  labs(title = "Log10 abundance of Actinobaculum (unique-to-semen) by site",
       subtitle = "Actinobaculum is more likely to be found in boar at Willingham",
       y = "log10 Relative Abundance")
ggsave("Actino_site.pdf", width = 8, height = 6, units = "in")
# Wow
kruskal.test(Actinobaculum_abund$log10_abund~Actinobaculum_abund$depth_binned)
kruskal.test(Actinobaculum_abund$log10_abund~Actinobaculum_abund$Line)
wilcox.test(Actinobaculum_abund$log10_abund~Actinobaculum_abund$Site)
# Quick comparison by Site
ggplot(Actinobaculum_abund, aes(x = Line, y = log10_abund, fill = Line)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  scale_y_continuous(labels = scales::math_format(10^.x)) +
  theme_classic() +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  scale_fill_igv() +
  labs(title = "Log10 Abundance of Actinobaculum (unique-to-semen) by Line")
# Not significant - not worth plotting.
# curious about read depth Fastidiosipila
# Extract Bacilli abundance
Fastidiosipila_abund <- as.data.frame(log_abund_matrix["Fastidiosipila", ])
Fastidiosipila_abund$Sample <- rownames(Fastidiosipila_abund)
colnames(Fastidiosipila_abund)[1] <- "log10_abund"
# Add metadata
Fastidiosipila_abund <- left_join(Fastidiosipila_abund, Metadata_semen, by = c("Sample" = "sampleid"))

# Quick comparison by Read depth
ggplot(Fastidiosipila_abund, aes(x = depth_binned, y = log10_abund, fill = depth_binned)) +
  geom_boxplot(width = 0.3, show.legend = F) +
  geom_jitter(width = 0.15, show.legend = F) +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  scale_y_continuous(labels = scales::math_format(10^.x)) +
  theme_classic(base_size = 10) +
  scale_fill_igv() +
  theme(text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")) +
  labs(title = "Log10 Abundance of Fastidiosipila (unique-to-semen) by Read Depth",
       y = "log10 Relative Abundance")
# Confusing as absent from med depth group - but some presence in low depth group
kruskal.test(Fastidiosipila_abund$log10_abund~Fastidiosipila_abund$depth_binned)
fastmodel <- lm(log10_abund~read_depth, data = Fastidiosipila_abund)
summary(fastmodel)
kruskal.test(Fastidiosipila_abund$log10_abund~Fastidiosipila_abund$Line)
wilcox.test(Fastidiosipila_abund$log10_abund~Fastidiosipila_abund$Site)
cor.test(Fastidiosipila_abund$log10_abund, Fastidiosipila_abund$read_depth, method = "spearman")

# Create cladograms
library(ggtree)
library(treeio)
library(ape)
library(phytools)
library(ggsci)
# Extract tree and taxonomy - the tree is already within the phyloseq object
# Qiime trees tend to be in Newick format
tree <- phy_tree(genus_glom)
tax <- as.data.frame(tax_table(genus_glom))
tax$label <- rownames(tax) #Ensure tip labels are a column called 'label'

# Create a presence/absence matrix for each niche
otu_mat <- otu_table(genus_glom)
otu_mat[otu_mat > 0] <- 1
otu_mat <- as.data.frame(t(otu_mat)) # Tranpose to taxa as rows
otu_mat$sampletype <- sample_data(genus_glom)$sampletype
# Summarize presence/absence by group
group_pa <- otu_mat |>
  group_by(sampletype) |>
  summarise(across(where(is.numeric), max)) |>
  column_to_rownames("sampletype") |>
  t() |>
  as.data.frame()

group_pa$label <- rownames(group_pa)

# Merge presence/absence with taxonomy
merged_df <- merge(group_pa, tax[,
                                 c("label",
                                   "Genus",
                                   "Family",
                                   "Order",
                                   "Class",
                                   "Phylum")],
                   by = "label",
                   all.x = TRUE)

# Assume tax is your taxonomy table as a data.frame
get_best_label <- function(row) {
  genus <- row["Genus"]
  family <- row["Family"]
  order <- row["Order"]
  class <- row["Class"]
  phylum <- row["Phylum"]
  if (!is.na(genus) && genus != "" && genus != "uncultured") {
    return(genus)
  } else if (!is.na(family) && family != "" && family != "uncultured") {
    return(paste0("Family: ", family))
  } else if (!is.na(order) && order != "" && order != "uncultured") {
    return(paste0("Order: ", order))
  } else if (!is.na(class) && class != "" && class != "uncultured") {
    return(paste0("Class: ", class))
  } else if (!is.na(phylum) && phylum != "" && phylum != "uncultured") {
    return(paste0("Phylum: ", phylum))
  } else {
    return("Unclassified")
  }
}
# Apply label function
merged_df$Best_Taxonomy <- apply(
  merged_df[, c("Genus", "Family", "Order", "Class", "Phylum")],
  1, get_best_label
)
# Assign niche categories
merged_df[, c("semen", "gloveswab", "extender", "boarswab")][
  is.na(merged_df[, c("semen", "gloveswab", "extender", "boarswab")])
] <- 0
# Combine niche presence into one variable
merged_df$niche <- apply(
  merged_df[, c("semen", "gloveswab", "extender", "boarswab")], 1,
  function(x) {
    present <- names(which(x == 1))
    if (length(present) == 0) {
      return("Absent")
    } else {
      return(paste(present, collapse = "_"))
    }
  }
)
# Try plotting semen only.
# Filter to only taxa unique to semen
semen_only_df <- merged_df[merged_df$niche == "semen", ]
# Subset tree to just these tips
semen_tree <- keep.tip(tree, semen_only_df$label)
# Generate best label (Genus > Family > Order > Class)
semen_only_df$Best_Taxonomy <- apply(
  semen_only_df[, c("Genus","Family","Order","Class","Phylum")],
  1, get_best_label
)
# ─────────────────────────────────────────────
# 2.  Create a 'highlight' flag for the two key genera
# ─────────────────────────────────────────────
key_genera <- c("Actinobaculum","Fastidiosipila")
semen_only_df <- semen_only_df %>%
  mutate(highlight = ifelse(Genus %in% key_genera, "key", "other"),
         shape     = ifelse(highlight == "key", 24, 21))   # 24 = filled triangle, 21 = circle
# Shorten long taxonomic names
semen_only_df$Best_Taxonomy <- sub("^Family: ", "", semen_only_df$Best_Taxonomy)
# ─────────────────────────────────────────────
# 3.  Plot: colour by Phylum, red outline for key genera
# ─────────────────────────────────────────────
# Build the tree object with tip positions
p_base <- ggtree(semen_tree, layout = "rectangular") %<+% semen_only_df

# Extract the plotted data (has x, y, everything)
tree_df <- p_base$data
# (after building tree_df with the numeric 'shape' column)
p <- p_base +
  geom_tippoint(aes(fill = Phylum, shape = shape),
                size = 2.5, colour = "black", stroke = 0.3) +
  scale_shape_identity(guide = "none") +                # interprets numeric 21 / 24
  scale_fill_igv(name = "Phylum",                               # coloured legend
                 guide = guide_legend(override.aes = list(shape = 21,
                                                          size  = 4,
                                                          colour = "black"))) +
  geom_text_repel(
    data = tree_df %>%
      filter(label %in% semen_only_df$label[semen_only_df$Genus %in% key_genera]),
    aes(x = x, y = y, label = Genus),
    size = 3, nudge_x = 0.4, min.segment.length = 0,
    inherit.aes = FALSE
  ) +
  theme_tree2() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7),
        text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")) +
  labs(title = "Phylogram of taxa unique to semen",
       subtitle = "The two most prevalent taxa in this group have been highlighted")
# after you’ve finished all other layers …
max_x <- max(tree_df$x)

p <- p +
  geom_treescale(
    x      = 0,          # x-coordinate where the bar starts
    y      = -1.5,         # y just below the lowest tip (negative puts it under the tree)
    width  = 0.1,        # branch-length represented by the bar (0.1 = 0.1 substitutions/site)
    offset = 0,          # extra gap between bar and tree (0 = flush)
    fontsize = 3
  ) +
  coord_cartesian(xlim = c(0, max_x), clip = "off") +
  theme(plot.margin = margin(1, 1, 5, 1, "pt"))  # extra bottom space for the bar
p
!is.null(semen_tree$edge.length) && all(!is.na(semen_tree$edge.length))
ggsave(plot = p, "semen_unique_clade.pdf", width = 9, height = 7)
ggsave(plot = p, "semen_unique_clade.png", width = 9, height = 7, dpi = 300)
# ─────────────────────────────────────────
# 1.  Subset metadata & tree
# ─────────────────────────────────────────
shared_df  <- merged_df %>% filter(niche != "semen")           # any tip found in ≥2 sample types
shared_tree <- ape::keep.tip(tree, shared_df$label)
absent_order <- setdiff(unique(shared_df$Order), unique(semen_only_df$Order))
print(absent_order)
# ─────────────────────────────────────────
# 1.  For each absent phylum, collect tip labels and MRCA node
# ─────────────────────────────────────────
hilite_nodes <- list()

for (order in absent_order) {
  order_tips <- shared_df$label[shared_df$Order == order]
  if (length(order_tips) > 1) {                       # need ≥2 tips for MRCA
    node  <- MRCA(shared_tree, order_tips)
    hilite_nodes[[order]] <- node
  }
}

# ─────────────────────────────────────────
# 2.  Flag Top-20 genera categories
#     (edit these vectors as you like)
# ─────────────────────────────────────────
likely_semen <- c("Dolosicoccus", "Lactobacillus", "Family: Carnobacteriaceae","Family: Aerococcaceae")
likely_env   <- c("Aliicoccus", "Streptococcus", "Staphylococcus", "Escherichia-Shigella", "Clostridium_sensu_stricto_1", "Acinetobacter")  # example
unknown <- c("Peptoniphilus",
             "Terrisporobacter", 
             "Corynebacterium", 
             "Psychrobacter",
             "W5053",
             "HT002",
             "Turicibacter")
shared_df <- shared_df |>
  mutate(origin = case_when(
    Best_Taxonomy %in% likely_semen ~ "Semen-origin",
    Best_Taxonomy %in% likely_env   ~ "Environmental",
    Best_Taxonomy %in% unknown ~ "Ambiguous",
    TRUE                             ~ "Unknown"),
    shape  = case_when(
      origin == "Semen-origin" ~ 24,   # filled triangle
      origin == "Environmental" ~ 22,  # filled square
      origin == "Ambiguous" ~ 21,
      TRUE                      ~ 20))  # circle
# 1. Make `origin` a factor in the order you want
shared_df$origin <- factor(shared_df$origin,
                           levels = c("Semen-origin", "Environmental", "Ambiguous", "Unknown"))

# 2. Shape palette keyed by the same origin names
shape_pal <- c("Semen-origin"  = 24,  # filled triangle
               "Environmental" = 22,  # filled square
               "Ambiguous"     = 21,  # filled circle
               "Unknown"       = 20)  # small dot
# 3. Build the tree
p_base <- ggtree(shared_tree, layout = "rectangular") %<+% shared_df
tree_df <- p_base$data

# 4. Plot
p_shared <- p_base +
  geom_tippoint(aes(fill = Phylum, shape = origin),   # map shape to *origin*, not to numeric
                size = 2.5, colour = "black", stroke = 0.35) +
  scale_shape_manual(name = "Origin", values = shape_pal) +
  scale_fill_igv(name = "Phylum",
                 guide = guide_legend(override.aes = list(shape = 21,
                                                          size  = 4,
                                                          colour = "black"))) +
  geom_tiplab(data   = tree_df %>% filter(origin != "Unknown"),
              aes(label = gsub("^Family: ", "", Best_Taxonomy)),
              align   = TRUE,
              linetype = "dotted",  # make leader lines subtle
              linesize = 0.25,
              offset = 0.2,
              size     = 2.8,
              hjust    = 0) +       # label text left-justified
  theme_tree2() +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7),
        theme(plot.margin = margin(1, 1, 3, 1, "pt")),
        title = element_text(size = 14, face = "bold"))
for (order in names(hilite_nodes)) {
  p_shared <- p_shared +
    geom_hilight(node = hilite_nodes[[order]],
                 fill = "grey50", alpha = 0.10, extend = 0.3)
}

for (order in names(hilite_nodes)) {
    p_shared <- p_shared +
      geom_cladelabel(
        node        = hilite_nodes[[order]],
        label       = order,
        align       = TRUE,
        offset      = 0.01,      # ← shorter bar
        offset.text = 0.025,      # ← text nudged same amount
        color       = "grey50",
        barsize     = 0.2,
        fontsize    = 3
      )
}
max_x <- max(tree_df$x)

p_shared <- p_shared +
  geom_treescale(
    x      = 0,          # x-coordinate where the bar starts
    y      = -1.5,         # y just below the lowest tip (negative puts it under the tree)
    width  = 0.1,        # branch-length represented by the bar (0.1 = 0.1 substitutions/site)
    offset = 0,          # extra gap between bar and tree (0 = flush)
    fontsize = 3
  ) +
  coord_cartesian(
    xlim  = c(0, max_x + 0.2),   # add generous space for text
    clip  = "off") + # allow grobs outside panel area
  labs(title = "Phylogram of shared taxa")

p_shared 
# right margin extra
ggsave("shared_taxa_cladogram_labeled.png", p_shared,
       width = 11, height = 9, dpi = 300)
ggsave("shared_taxa_cladogram_labeled.pdf", p_shared,
       width = 11, height = 9)

# Create a phyloseq object with the 4 comparative samples plus mocks
ps_amplicon <- subset_samples(ps_filtered,
                               (sample_names(ps_filtered) %in% c("42A0118K",
                                                                 "5A0295K",
                                                                 "2A0297K",
                                                                 "41A0888K",
                                                                 "39A9968H",
                                                                 "50ADNAMCII")))
# Save objects that we'll need for the shotgun and RNA scripts.
saveRDS(ps_amplicon, "shotgun/ps_amplicon.rds")
