setwd("~/Desktop/InducibleArchitectures/raw_data/")

library("ggplot2")
library("dplyr")
library("tidyr")
require(cowplot)

options(stringsAsFactors = F)

# GET RAW BARCODE READS AND READ IN FILES
filelist = list.files(pattern = 'fLP3*') # change prefix for all barcode counts
for(i in filelist) {
  x <- read.table(i, col.names=c(i, 'barcode'), header = F)
  assign(i,x) 
}

# function that creates dataframe with 3 columns. 1 with raw reads, 1 normalized, 1 barcode.
make_normalized_df <- function(df) {
  norm = "normalized_"
  new_name = paste(norm, deparse(substitute(df)), sep = '')
  df <- df %>%
    mutate(temp_name = 1000000*df[[1]]/sum(df[[1]]))
  names(df)[names(df) == 'temp_name'] = new_name
  return(df)
}

# Make all normalized dataframes for each sample
fLP3_LacZ_Induced_DNA_1.txt <- make_normalized_df(fLP3_LacZ_Induced_DNA_1.txt)
fLP3_LacZ_Induced_DNA_2.txt <- make_normalized_df(fLP3_LacZ_Induced_DNA_2.txt)
fLP3_LacZ_Induced_RNA_1.txt <- make_normalized_df(fLP3_LacZ_Induced_RNA_1.txt)
fLP3_LacZ_Induced_RNA_2.txt <- make_normalized_df(fLP3_LacZ_Induced_RNA_2.txt)

fLP3_LacZ_UnInduced_DNA_1.txt <- make_normalized_df(fLP3_LacZ_UnInduced_DNA_1.txt)
fLP3_LacZ_UnInduced_DNA_2.txt <- make_normalized_df(fLP3_LacZ_UnInduced_DNA_2.txt)
fLP3_LacZ_UnInduced_RNA_1.txt <- make_normalized_df(fLP3_LacZ_UnInduced_RNA_1.txt)
fLP3_LacZ_UnInduced_RNA_2.txt <- make_normalized_df(fLP3_LacZ_UnInduced_RNA_2.txt)

# create combined dataframe called flp3_loop
fLP3_loop <- full_join(fLP3_LacZ_Induced_DNA_1.txt, fLP3_LacZ_Induced_DNA_2.txt, by='barcode') %>% # change .txt files to files you've read in
  full_join(., fLP3_LacZ_Induced_RNA_1.txt, by='barcode') %>%
  full_join(., fLP3_LacZ_Induced_RNA_2.txt, by='barcode') %>%
  full_join(., fLP3_LacZ_UnInduced_DNA_1.txt, by='barcode') %>%
  full_join(., fLP3_LacZ_UnInduced_DNA_2.txt, by='barcode') %>%
  full_join(., fLP3_LacZ_UnInduced_RNA_1.txt, by='barcode') %>%
  full_join(., fLP3_LacZ_UnInduced_RNA_2.txt, by='barcode')

names(fLP3_loop) = sub(".txt","", names(fLP3_loop)) #rename all colummns that were named after text file

rm(list = c(filelist))
rm(x)

# Read in all barcodes from barcode stats file
barcode_stats_loop <- read.table("LacZ_mapping_barcode_statistics.txt", header = T)
mapped_barcodes <- barcode_stats_loop[!is.na(barcode_stats_loop$most_common),] # removed unmapped
Compare_barcode_Reps <- left_join(mapped_barcodes, fLP3_loop, by ='barcode') # left join with sequencing data
Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0 #Make all NA values = 0
temp <- Compare_barcode_Reps # save it in temp dataframe 

#Jess worked on this section --------------------------------------

var <- temp %>%
  select(-num_unique, -num_reads, -num_reads_most_common)

induced_bc_1 <- var %>%
  select(most_common, barcode, fLP3_LacZ_Induced_DNA_1, normalized_fLP3_LacZ_Induced_DNA_1, fLP3_LacZ_Induced_RNA_1, normalized_fLP3_LacZ_Induced_RNA_1) %>%
  filter(fLP3_LacZ_Induced_DNA_1 > 2) %>%
  group_by(most_common) %>%
  mutate(median_RNA_exp_Induced_1 = median(normalized_fLP3_LacZ_Induced_RNA_1/normalized_fLP3_LacZ_Induced_DNA_1)) %>%
  mutate(sum_RNA_exp_Induced_1 = sum(normalized_fLP3_LacZ_Induced_RNA_1)/sum(normalized_fLP3_LacZ_Induced_DNA_1)) %>%
  ungroup() %>%
  select(most_common, median_RNA_exp_Induced_1, sum_RNA_exp_Induced_1) %>%
  distinct()

induced_bc_2 <- var %>%
  select(most_common, barcode, fLP3_LacZ_Induced_DNA_2, normalized_fLP3_LacZ_Induced_DNA_2, fLP3_LacZ_Induced_RNA_2, normalized_fLP3_LacZ_Induced_RNA_2) %>%
  filter(fLP3_LacZ_Induced_DNA_2 > 2) %>%
  group_by(most_common) %>%
  mutate(sum_RNA_exp_Induced_2 = sum(normalized_fLP3_LacZ_Induced_RNA_2)/sum(normalized_fLP3_LacZ_Induced_DNA_2)) %>%
  ungroup() %>%
  select(most_common, median_RNA_exp_Induced_2, sum_RNA_exp_Induced_2) %>%
  distinct()

uninduced_bc_1 <- var %>%
  select(most_common, barcode, fLP3_LacZ_UnInduced_DNA_1, normalized_fLP3_LacZ_UnInduced_DNA_1, fLP3_LacZ_UnInduced_RNA_1, normalized_fLP3_LacZ_UnInduced_RNA_1) %>%
  filter(fLP3_LacZ_UnInduced_DNA_1 > 2) %>%
  group_by(most_common) %>%
  mutate(median_RNA_exp_UnInduced_1 = median(normalized_fLP3_LacZ_UnInduced_RNA_1/normalized_fLP3_LacZ_UnInduced_DNA_1)) %>%
  mutate(sum_RNA_exp_UnInduced_1 = sum(normalized_fLP3_LacZ_UnInduced_RNA_1)/sum(normalized_fLP3_LacZ_UnInduced_DNA_1)) %>%
  ungroup() %>%
  select(most_common, median_RNA_exp_UnInduced_1, sum_RNA_exp_UnInduced_1) %>%
  distinct()


uninduced_bc_2 <- var %>%
  select(most_common, barcode, fLP3_LacZ_UnInduced_DNA_2, normalized_fLP3_LacZ_UnInduced_DNA_2, fLP3_LacZ_UnInduced_RNA_2, normalized_fLP3_LacZ_UnInduced_RNA_2) %>%
  filter(fLP3_LacZ_UnInduced_DNA_2 > 2) %>%
  group_by(most_common) %>%
  mutate(median_RNA_exp_UnInduced_2 = median(normalized_fLP3_LacZ_UnInduced_RNA_2/normalized_fLP3_LacZ_UnInduced_DNA_2)) %>%
  mutate(sum_RNA_exp_UnInduced_2 = sum(normalized_fLP3_LacZ_UnInduced_RNA_2)/sum(normalized_fLP3_LacZ_UnInduced_DNA_2)) %>%
  ungroup() %>%
  select(most_common, median_RNA_exp_UnInduced_2, sum_RNA_exp_UnInduced_2) %>%
  distinct()

# combine DF's at variant level
variant_stats <- read.table("LacZ_mapping_variant_statistics.txt", header = T, fill = T)

induce_exp <- variant_stats %>% select(name, variant) %>%
  left_join(., induced_bc_1, by = c("variant" = "most_common")) %>%
  left_join(., induced_bc_2, by = c("variant" = "most_common")) %>%
  left_join(., uninduced_bc_1, by = c("variant" = "most_common")) %>%
  left_join(., uninduced_bc_2, by = c("variant" = "most_common")) %>%
  na.omit(.)

# store all the negative controls
test <- subset(induce_exp, grepl("neg_control", induce_exp$name)) 

# median of median neg_control expression
med_neg_exp_ind_1 <- median(induce_exp$median_RNA_exp_Induced_1)
med_neg_exp_ind_2 <- median(test$median_RNA_exp_Induced_2)
med_neg_exp_un_1 <- median(test$median_RNA_exp_UnInduced_1)
med_neg_exp_un_2 <- median(test$median_RNA_exp_UnInduced_2)

# median of sum neg_control expression
sum_neg_exp_ind_1 <- median(induce_exp$sum_RNA_exp_Induced_1)
sum_neg_exp_ind_2 <- median(test$sum_RNA_exp_Induced_2)
sum_neg_exp_un_1 <- median(test$sum_RNA_exp_UnInduced_1)
sum_neg_exp_un_2 <- median(test$sum_RNA_exp_UnInduced_2)

# normalize to median or sum neg_control expression
induce_exp <- induce_exp %>%
  mutate(median_normalized_RNA_exp_Induced_1 = median_RNA_exp_Induced_1/med_neg_exp_ind_1,
         median_normalized_RNA_exp_Induced_2 = median_RNA_exp_Induced_2/med_neg_exp_ind_2,
         median_normalized_RNA_exp_UnInduced_1 = median_RNA_exp_UnInduced_1/med_neg_exp_un_1,
         median_normalized_RNA_exp_UnInduced_2 = median_RNA_exp_UnInduced_2/med_neg_exp_un_2,
         median_normalized_RNA_exp_Induced_12 = (median_normalized_RNA_exp_Induced_1 + median_normalized_RNA_exp_Induced_2)/2,
         median_normalized_RNA_exp_UnInduced_12 = (median_normalized_RNA_exp_UnInduced_1 + median_normalized_RNA_exp_UnInduced_2)/2,
         median_ratio = median_normalized_RNA_exp_Induced_12/median_normalized_RNA_exp_UnInduced_12) %>%
  mutate(normalized_RNA_exp_Induced_1 = sum_RNA_exp_Induced_1/sum_neg_exp_ind_1,
         normalized_RNA_exp_Induced_2 = sum_RNA_exp_Induced_2/sum_neg_exp_ind_2,
         normalized_RNA_exp_UnInduced_1 = sum_RNA_exp_UnInduced_1/sum_neg_exp_un_1,
         normalized_RNA_exp_UnInduced_2 = sum_RNA_exp_UnInduced_2/sum_neg_exp_un_2,
         normalized_RNA_exp_Induced_12 = (normalized_RNA_exp_Induced_1 + normalized_RNA_exp_Induced_2)/2,
         normalized_RNA_exp_UnInduced_12 = (normalized_RNA_exp_UnInduced_1 + normalized_RNA_exp_UnInduced_2)/2,
         ratio = normalized_RNA_exp_Induced_12/normalized_RNA_exp_UnInduced_12) %>%
  filter(normalized_RNA_exp_UnInduced_12 > 0) %>%
  select(name, variant, normalized_RNA_exp_Induced_1, normalized_RNA_exp_Induced_2, normalized_RNA_exp_UnInduced_1,
         normalized_RNA_exp_UnInduced_2, normalized_RNA_exp_Induced_12, normalized_RNA_exp_UnInduced_12, ratio)

library(GGally)
# correlation plots
comp <- read.table('~/Desktop/newFlowComparison.txt', header = TRUE)
comp <- comp %>% select(-Name)
ggpairs(comp, columnLabels = c("Flow", "Median", "Sum"), lower = list(
  continuous = "smooth"))


cor(induce_exp$median_ratio, induce_exp$sum_ratio) #Pearson correlation
induce_exp %>%
  ggplot(aes(x = log10(median_ratio), y = log10(sum_ratio))) +
  geom_point(alpha = 0.75, fill = 'black')


#You could try following this format for all variants where each sample is
#processed separately. Maintaining the DNA samples with their corresponding RNA
#replicates (I thought these were taken from the same flasks anyway). Then you 
#would combine all dfs at the variant level, taking the average of 
#sum or median_ratio_norm between replicates, then determine the induction ratio 
#with this average_ind/average_unind. If you do this you would also get slightly
#different values for the negative controls that you normalize to so the values 
#I calculated here for S1 may end up being slightly different.

#End Jess section---------------------------------------------------------------


# FUNCTION FOR CREATING DATAFRAME MADE WITH SUMMED EXP, num_barcodes > 2 filter
sum_var <- function(df) {
  uninduced <- df %>%
    select(most_common, normalized_fLP3_LacZ_UnInduced_DNA_1, normalized_fLP3_LacZ_UnInduced_DNA_2, normalized_fLP3_LacZ_UnInduced_RNA_1,
           normalized_fLP3_LacZ_UnInduced_RNA_2) %>%
    # filter(normalized_fLP3_LacZ_UnInduced_DNA_1 > 0 | normalized_fLP3_LacZ_UnInduced_DNA_2 > 0) %>%
    group_by(most_common) %>%
    mutate(num_barcodes_unind = n()) %>%
    filter(num_barcodes_unind > 2) %>%
    mutate(DNA_UnInduced_avg = mean(sum(normalized_fLP3_LacZ_UnInduced_DNA_1), sum(normalized_fLP3_LacZ_UnInduced_DNA_2))) %>%
    mutate(RNA_exp_UnInduced_1 = sum(normalized_fLP3_LacZ_UnInduced_RNA_1)/(DNA_UnInduced_avg),
           RNA_exp_UnInduced_2 = sum(normalized_fLP3_LacZ_UnInduced_RNA_2)/(DNA_UnInduced_avg)) %>%
    ungroup() %>%
    select(-normalized_fLP3_LacZ_UnInduced_DNA_1, -normalized_fLP3_LacZ_UnInduced_DNA_2, -normalized_fLP3_LacZ_UnInduced_RNA_1, -normalized_fLP3_LacZ_UnInduced_RNA_2) %>%
    distinct()
  
  induced <- df %>%
    select(most_common, normalized_fLP3_LacZ_Induced_DNA_1, normalized_fLP3_LacZ_Induced_DNA_2, normalized_fLP3_LacZ_Induced_RNA_1,
           normalized_fLP3_LacZ_Induced_RNA_2) %>%
    # filter(normalized_fLP3_LacZ_Induced_DNA_1 > 0  | normalized_fLP3_LacZ_Induced_DNA_2 > 0) %>%
    group_by(most_common) %>%
    mutate(num_barcodes_ind = n()) %>%
    filter(num_barcodes_ind > 2) %>%
    mutate(DNA_Induced_avg = mean(sum(normalized_fLP3_LacZ_Induced_DNA_1), sum(normalized_fLP3_LacZ_Induced_DNA_2))) %>%
    mutate(RNA_exp_Induced_1 = sum(normalized_fLP3_LacZ_Induced_RNA_1)/(DNA_Induced_avg),
           RNA_exp_Induced_2 = sum(normalized_fLP3_LacZ_Induced_RNA_2)/(DNA_Induced_avg)) %>%
    ungroup() %>%
    select(-normalized_fLP3_LacZ_Induced_DNA_1, -normalized_fLP3_LacZ_Induced_DNA_2, -normalized_fLP3_LacZ_Induced_RNA_1, -normalized_fLP3_LacZ_Induced_RNA_2) %>%
    distinct()
  inner_join(induced, uninduced, by = 'most_common')
}
mean_exp_master <- sum_var(temp) 

# read in variant file to get all names of promoters 
variant_stats <- read.table("LacZ_mapping_variant_statistics.txt", header = T, fill = T)

# create dataframe with all variant names and filter out expression > 0
induce_exp <- left_join(variant_stats, mean_exp_master, by = c("variant" = "most_common")) %>% 
  na.omit(.) %>% filter(RNA_exp_Induced_1 > 0 & RNA_exp_Induced_2 > 0 & RNA_exp_UnInduced_1 > 0 & RNA_exp_UnInduced_2 > 0) %>% 
  select("variant", "name", "RNA_exp_Induced_1", "RNA_exp_Induced_2", "RNA_exp_UnInduced_1", "RNA_exp_UnInduced_2", num_barcodes_ind, num_barcodes_unind)

# store all the negative controls
test <- subset(induce_exp, grepl("neg_control", induce_exp$name)) 

# median of neg_control expression
med_neg_exp_ind_1 <- median(test$RNA_exp_Induced_1)
med_neg_exp_ind_2 <- median(test$RNA_exp_Induced_2)
med_neg_exp_un_1 <- median(test$RNA_exp_UnInduced_1)
med_neg_exp_un_2 <- median(test$RNA_exp_UnInduced_2)

# for the negative controls
test <- test %>%
  mutate(normalized_RNA_exp_Induced_1 = RNA_exp_Induced_1/med_neg_exp_ind_1,
         normalized_RNA_exp_Induced_2 = RNA_exp_Induced_2/med_neg_exp_ind_2,
         normalized_RNA_exp_UnInduced_1 = RNA_exp_UnInduced_1/med_neg_exp_un_1,
         normalized_RNA_exp_UnInduced_2 = RNA_exp_UnInduced_2/med_neg_exp_un_2,
         normalized_RNA_exp_Induced_12 = (normalized_RNA_exp_Induced_1 + normalized_RNA_exp_Induced_2)/2,
         normalized_RNA_exp_UnInduced_12 = (normalized_RNA_exp_UnInduced_1 + normalized_RNA_exp_UnInduced_2)/2,
         ratio = normalized_RNA_exp_Induced_12/normalized_RNA_exp_UnInduced_12)

# expression for each sample is normalized first by dividing by the median of the negative control expression in that c
# orresponding sample. These normalized expressions are then averaged across the two replicates for Induced and UnInduced. 
# Ratio is determined by dividing normalized Induced expression by normalized UnInduced expression. 
induce_exp <- induce_exp %>%
  mutate(normalized_RNA_exp_Induced_1 = RNA_exp_Induced_1/med_neg_exp_ind_1,
         normalized_RNA_exp_Induced_2 = RNA_exp_Induced_2/med_neg_exp_ind_2,
         normalized_RNA_exp_UnInduced_1 = RNA_exp_UnInduced_1/med_neg_exp_un_1,
         normalized_RNA_exp_UnInduced_2 = RNA_exp_UnInduced_2/med_neg_exp_un_2,
         normalized_RNA_exp_Induced_12 = (normalized_RNA_exp_Induced_1 + normalized_RNA_exp_Induced_2)/2,
         normalized_RNA_exp_UnInduced_12 = (normalized_RNA_exp_UnInduced_1 + normalized_RNA_exp_UnInduced_2)/2,
         ratio = normalized_RNA_exp_Induced_12/normalized_RNA_exp_UnInduced_12)

induce_exp$name <- gsub("shift", "", induce_exp$name)
induce_exp$name <- gsub(">LACI-", "", induce_exp$name)
induce_exp$name <- gsub("-flipped_rc", "", induce_exp$name)
induce_exp$name <- gsub("LOOP", "", induce_exp$name)
induce_exp$name <- gsub("bp", "", induce_exp$name)

backup <- induce_exp
write.table(backup, "../processed_data/induce_exp(updated).txt", quote = F, row.names = F)
#induce_exp <- read.table("./induce_exp.txt", header = T)
