# 20/01/2025##
##By: Octavio Zambada##

# This script processes multiple input files to compare and merge data, generating a final output file named
# "HormoneOrthologs.txt". This file compiles information on three different phytohormones in tomato and their 
# orthologs in both wild and domesticated species.
# The script applies various functions to:
# 1. Identify and compare orthologs across species.
# 2. Integrate relevant phytohormone-related data.
# 3. Output a structured, tab-delimited text file containing the compiled results.
# 4. Generate output files to enhance network visualizations.

## Set working directory
setwd("C:/Users/octav/OneDrive/Escritorio/TomatoOrtho/analisis2/")

# Load MYC2 interactors file
MYC <- read.table("C:/Users/octav/OneDrive/Escritorio/TomatoOrtho/MYCinteractorsEdit.txt",
                  quote="\"", comment.char="")

# Load orthologs file
orthologs <- read.delim("C:/Users/octav/OneDrive/Escritorio/TomatoOrtho/orthologs_sinPunto.tsv")

# Load dplyr library for data manipulation
library(dplyr)

#### Join both datasets
# Perform an inner join between MYC data and orthologs
result <- inner_join(MYC, orthologs, by = c("V1" = "b"))
result2 <- inner_join(MYC, orthologs, by = c("V1" = "a"))

# Rename column "b" to "a" in the second result
colnames(result2)[colnames(result2) == "b"] <- "a"

# Combine both results into one dataframe
joinMYC <- rbind(result, result2)

# Remove rows with missing values
joinMYC <- na.omit(joinMYC)

# Check how many unique values exist in column "a"
length(unique(joinMYC$a))

# Reload dplyr to ensure the library is available
library(dplyr)

# Rename columns for clarity
colnames(joinMYC) <- c("Gen_a", "Gen_b", "Species_b", "Species_a",
                       "OG", "Normalized_bitScore")

# Create a frequency table for tomato genes
tableMYC <- as.data.frame(table(joinMYC$Gen_a))

# Replace "tomato" in Species_b with the corresponding value in Species_a
joinMYC <- joinMYC %>%
  mutate(Species_b = ifelse(Species_b == "tomato" & !is.na(Species_b), Species_a, Species_b))

# Change all values in Species_a to "tomato", but keep NAs in that column
joinMYC <- joinMYC %>%
  mutate(Species_a = ifelse(!is.na(Species_a), "tomato", Species_a))

# Add hormone type as "Jasmonic Acid"
joinMYC$Hormone <- "Jasmonic Acid"

# Rename columns again to include the hormone type
colnames(joinMYC) <- c("Gen_a", "Gen_b", "Species_b", "Species_a",
                       "OG", "Normalized_bitScore", "Hormone")

# Count how many times each tomato gene appears in other species
duplicated_genesMYC <- joinMYC %>%
  group_by(Gen_a, Species_b) %>%  # Group by gene and species
  summarise(count = n(), .groups = "drop") %>%  # Count occurrences
  filter(count > 1)

# Load MYC2 network data
MYC2net <- read.delim("C:/Users/octav/OneDrive/Escritorio/TomatoOrtho/MYC2_net.txt.sif", 
                      header=FALSE)

# Clean the MYC2 network data by removing suffixes
MYC2net$V1 <- gsub("\\.1$", "", MYC2net$V1)
MYC2net$V3 <- gsub("\\.1$", "", MYC2net$V3)
MYC2net$V1 <- gsub("\\.2$", "", MYC2net$V1)
MYC2net$V3 <- gsub("\\.2$", "", MYC2net$V3)

# Merge the datasets, keeping all rows from MYC2net
merged_dataMYC2 <- merge(MYC2net, tableMYC, by.x = "V3", by.y = "Var1", all.x = TRUE)

# Replace NA in the Freq column with 0
merged_dataMYC2$Freq[is.na(merged_dataMYC2$Freq)] <- 0

# Rename the "Freq" column to "Orthologs"
colnames(merged_dataMYC2)[colnames(merged_dataMYC2) == "Freq"] <- "Ortologs"

# Extract the necessary columns
MYC2numOrto <- merged_dataMYC2[, c(1, 4)]

# Rename the columns for clarity
colnames(MYC2numOrto) <- c("Gen", "N.Orthologs")

# Write the MYC2 orthologs data to a file
write.table(MYC2numOrto, file = "MYC2numOrtho.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Write the join data to a file
write.table(joinMYC, file = "MYC2OrthoSpecies.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Write the merged MYC2 network data to a file
write.table(merged_dataMYC2, file = "MYC2_net.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Write the duplicated genes data to a file
write.table(duplicated_genesMYC , file = "duplicated_genesSpeciesMYC2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(unique(duplicated_genesMYC$Gen_a)), file = "duplicatedMYC2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Brassinosteroids

# Read the BrassinoMip dataset
BrassinoMip <- read.table("C:/Users/octav/OneDrive/Escritorio/TomatoOrtho/BrassinoMipedit.txt",
                          quote="\"", comment.char="")

#### Join both datasets
# Perform a left join between BrassinoMip and orthologs
resultMip <- left_join(BrassinoMip, orthologs, by = c("V1" = "b"))
resultMip2 <- left_join(BrassinoMip, orthologs, by = c("V1" = "a"))

# Rename column 'b' to 'a' in the second result
colnames(resultMip2)[colnames(resultMip2) == "b"] <- "a"

# Combine the results into one dataset
joinMip <- rbind(resultMip, resultMip2)

# Check the number of unique genes
length(unique(joinMip$V1))

# Create a table of gene occurrences in the dataset
table(joinMip$Gen_a)

# Rename columns for clarity
colnames(joinMip) <- c("Gen_a", "Gen_b", "Species_b", "Species_a", "OG", "Normalized_bitScore")

# Write the joined dataset to a file
write.table(joinMip, file = "BrassinoMiPOrthoSpecies.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Replace "tomato" in Species_b with the corresponding value in Species_a
joinMip <- joinMip %>%
  mutate(Species_b = ifelse(Species_b == "tomato" & !is.na(Species_b), Species_a, Species_b))

# Change all values in Species_a to "tomato"
joinMip <- joinMip %>%
  mutate(Species_a = ifelse(!is.na(Species_a), "tomato", Species_a))

# Add hormone type as "Brassinosteroid"
joinMip$Hormone <- "Brassinosteroid"

# Create a table for BrassinoMip genes
tablemip <- as.data.frame(table(joinMip$Gen_a))

# Read the BrassinoMip network data
mipnet <- read.delim("C:/Users/octav/OneDrive/Escritorio/TomatoOrtho/BrassinoMip.txt.sif", 
                     header=FALSE)

# Clean the network data by removing suffixes
mipnet$V1 <- gsub("\\.1$", "", mipnet$V1)
mipnet$V3 <- gsub("\\.1$", "", mipnet$V3)
mipnet$V1 <- gsub("\\.2$", "", mipnet$V1)
mipnet$V3 <- gsub("\\.2$", "", mipnet$V3)

# Merge the network data with the table of gene occurrences, keeping all rows from mipnet
merged_datamip <- merge(mipnet, tablemip, by.x = "V3", by.y = "Var1", all.x = TRUE)

# Replace NA in the 'Freq' column with 0
merged_datamip$Freq[is.na(merged_datamip$Freq)] <- 0

# Rename the 'Freq' column to 'Ortologs'
colnames(merged_datamip)[colnames(merged_datamip) == "Freq"] <- "Ortologs"

# Extract the necessary columns
mipnumOrto <- merged_datamip[, c(1, 4)]

# Rename the columns for clarity
colnames(mipnumOrto) <- c("Gen", "N.Orthologs")

# Write the orthologs data to a file
write.table(mipnumOrto, file = "mipnumOrtho.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Write the join data to a file
write.table(joinMYC, file = "mipOrthoSpecies.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Write the merged network data to a file
write.table(merged_datamip, file = "mip_net.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Count how many times each tomato gene appears in other species
duplicated_genesMip <- joinMip %>%
  group_by(Gen_a, Species_b) %>%  # Group by gene and species
  summarise(count = n(), .groups = "drop") %>%  # Count occurrences
  filter(count > 1)

# Write the duplicated genes data to a file
write.table(duplicated_genesMip, file = "duplicated_genesSpeciesBrassinoMip.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(unique(duplicated_genesMip$Gen_a)), file = "duplicatedBrassinoMip.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


## AuxinMip

# Read the AuxinMip dataset
AuxinMip <- read.table("C:/Users/octav/OneDrive/Escritorio/TomatoOrtho/AuxinMipedit.txt",
                       quote="\"", comment.char="")

#### Join both datasets
# Perform a left join between AuxinMip and orthologs
resultAuxMip <- left_join(AuxinMip, orthologs, by = c("V1" = "b"))
resultAuxMip2 <- left_join(AuxinMip, orthologs, by = c("V1" = "a"))
resultAuxMip2 <- na.omit(resultAuxMip2)

# Rename column 'b' to 'a' in the second result
colnames(resultAuxMip2)[colnames(resultAuxMip2) == "b"] <- "a"

# Combine the results into one dataset
joinAuxMip <- rbind(resultAuxMip, resultAuxMip2)

# Check the number of unique genes
length(unique(joinAuxMip$V1))

# Create a table of gene occurrences in the dataset
table(joinAuxMip$V1)

# Create a table of species occurrences
table(joinAuxMip$species_a)
table(joinMYC$Species_a)

# Rename columns for clarity
colnames(joinAuxMip) <- c("Gen_a", "Gen_b", "Species_b", "Species_a", "OG", "Normalized_bitScore")

# Write the joined dataset to a file
write.table(joinAuxMip, file = "AuxMiPOrthoSpecies.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Replace "tomato" in Species_b with the corresponding value in Species_a
joinAuxMip <- joinAuxMip %>%
  mutate(Species_b = ifelse(Species_b == "tomato" & !is.na(Species_b), Species_a, Species_b))

# Change all values in Species_a to "tomato"
joinAuxMip <- joinAuxMip %>%
  mutate(Species_a = ifelse(!is.na(Species_a), "tomato", Species_a))

# Add hormone type as "Auxin"
joinAuxMip$Hormone <- "Auxin"

# Create a table for AuxinMip genes
tableAuxMip <- as.data.frame(table(joinAuxMip$Gen_a))

# Read the AuxinMip network data
AuxMipnet <- read.delim("C:/Users/octav/OneDrive/Escritorio/TomatoOrtho/AuxinMip.txt.sif", 
                        header=FALSE)

# Clean the network data by removing suffixes
AuxMipnet$V1 <- gsub("\\.1$", "", AuxMipnet$V1)
AuxMipnet$V3 <- gsub("\\.1$", "", AuxMipnet$V3)
AuxMipnet$V1 <- gsub("\\.2$", "", AuxMipnet$V1)
AuxMipnet$V3 <- gsub("\\.2$", "", AuxMipnet$V3)

# Merge the network data with the table of gene occurrences, keeping all rows from AuxMipnet
merged_dataAuxMip <- merge(AuxMipnet, tableAuxMip, by.x = "V3", by.y = "Var1", all.x = TRUE)

# Replace NA in the 'Freq' column with 0
merged_dataAuxMip$Freq[is.na(merged_dataAuxMip$Freq)] <- 0

# Rename the 'Freq' column to 'Ortologs'
colnames(merged_dataAuxMip)[colnames(merged_dataAuxMip) == "Freq"] <- "Ortologs"

# Extract the necessary columns
AuxMipnumOrto <- merged_dataAuxMip[, c(1, 4)]

# Rename the columns for clarity
colnames(AuxMipnumOrto) <- c("Gen", "N.Orthologs")

# Write the orthologs data to a file
write.table(AuxMipnumOrto, file = "AuxMipnumOrtho.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Write the merged network data to a file
write.table(merged_dataAuxMip, file = "AuxMip_net.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Count how many times each tomato gene appears in other species
duplicated_genesAuxMip <- joinAuxMip %>%
  group_by(Gen_a, Species_b) %>%  # Group by gene and species
  summarise(count = n(), .groups = "drop") %>%  # Count occurrences
  filter(count > 1)

# Write the duplicated genes data to a file
write.table(duplicated_genesAuxMip, file = "duplicated_genesSpeciesAuxMip.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(unique(duplicated_genesAuxMip$Gen_a)), file = "duplicatedAuxMip.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Combine the final data from all hormones
merged_data <- rbind(joinMYC, joinMip, joinAuxMip)

# Load libraries for visualization
library(ggplot2)
library(ggalluvial)
library(dplyr)

# Ensure the Ortholog_Status column is present and correct
merged_data <- merged_data %>%
  mutate(Ortholog_Status = ifelse(is.na(Gen_b), "No Ortholog", "Has Ortholog"),
         Species_with_Ortholog = ifelse(!is.na(Species_b), Species_b, "None"))

# Filter relevant columns
merged_filt <- merged_data[, c(1, 2, 8, 9, 7)]

# Identify genes with duplicated orthologs, excluding "None" and NA entries
merged_filt <- merged_filt %>%
  group_by(Gen_a, Species_with_Ortholog) %>%
  mutate(Duplicated_Ortholog = ifelse(n() > 1 & Species_with_Ortholog != "None" & !is.na(Species_with_Ortholog), TRUE, FALSE)) %>%
  ungroup()

# Write the final filtered data to a file
write.table(merged_filt, file = "HormoneOrthologs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Clean species names for visualization
merged_filt$Species_with_Ortholog <- gsub("\\.pep", "", merged_filt$Species_with_Ortholog)

# Convert Duplicated_Ortholog to factor (as before)
merged_filt$Duplicated_Ortholog <- factor(merged_filt$Duplicated_Ortholog, levels = c(TRUE, FALSE))

# Create an alluvial plot with 4 columns using geom_label
ggplot(data = merged_filt,
       aes(axis1 = Hormone, axis2 = Ortholog_Status, axis3 = Species_with_Ortholog, axis4 = Duplicated_Ortholog, y = 1)) +
  geom_alluvium(aes(fill = Hormone), width = 1/6) +
  geom_stratum(width = 1/6, fill = "white", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3,  # Use geom_label
             position = position_nudge(y = 0.05)) +  # Adjust vertical position with nudge
  scale_x_discrete(limits = c("Hormone", "Ortholog Status", "Species with Ortholog", "Duplicated Ortholog")) +
  labs(x = NULL,
       y = "Number of Genes") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Auxin" = "hotpink2", "Brassinosteroid" = "lightblue2", "Jasmonic Acid" = "palegreen2")) # Pastel colors

