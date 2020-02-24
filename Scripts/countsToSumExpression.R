setwd("../raw_data/alternative_architectures")

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
  mutate(num_barcodes_ind = n()) %>%
  filter(num_barcodes_ind > 2) %>%
  mutate(sum_RNA_exp_Induced_1 = sum(normalized_fLP3_LacZ_Induced_RNA_1)/sum(normalized_fLP3_LacZ_Induced_DNA_1)) %>%
  ungroup() %>%
  select(most_common, sum_RNA_exp_Induced_1, num_barcodes_ind) %>%
  distinct()

induced_bc_2 <- var %>%
  select(most_common, barcode, fLP3_LacZ_Induced_DNA_2, normalized_fLP3_LacZ_Induced_DNA_2, fLP3_LacZ_Induced_RNA_2, normalized_fLP3_LacZ_Induced_RNA_2) %>%
  filter(fLP3_LacZ_Induced_DNA_2 > 2) %>%
  group_by(most_common) %>%
  mutate(num_barcodes_ind = n()) %>%
  filter(num_barcodes_ind > 2) %>%
  mutate(sum_RNA_exp_Induced_2 = sum(normalized_fLP3_LacZ_Induced_RNA_2)/sum(normalized_fLP3_LacZ_Induced_DNA_2)) %>%
  ungroup() %>%
  select(most_common, sum_RNA_exp_Induced_2, num_barcodes_ind) %>%
  distinct()

uninduced_bc_1 <- var %>%
  select(most_common, barcode, fLP3_LacZ_UnInduced_DNA_1, normalized_fLP3_LacZ_UnInduced_DNA_1, fLP3_LacZ_UnInduced_RNA_1, normalized_fLP3_LacZ_UnInduced_RNA_1) %>%
  filter(fLP3_LacZ_UnInduced_DNA_1 > 2) %>%
  group_by(most_common) %>%
  mutate(num_barcodes_unind = n()) %>%
  filter(num_barcodes_unind > 2) %>%
  mutate(sum_RNA_exp_UnInduced_1 = sum(normalized_fLP3_LacZ_UnInduced_RNA_1)/sum(normalized_fLP3_LacZ_UnInduced_DNA_1)) %>%
  ungroup() %>%
  select(most_common, sum_RNA_exp_UnInduced_1, num_barcodes_unind) %>%
  distinct()


uninduced_bc_2 <- var %>%
  select(most_common, barcode, fLP3_LacZ_UnInduced_DNA_2, normalized_fLP3_LacZ_UnInduced_DNA_2, fLP3_LacZ_UnInduced_RNA_2, normalized_fLP3_LacZ_UnInduced_RNA_2) %>%
  filter(fLP3_LacZ_UnInduced_DNA_2 > 2) %>%
  group_by(most_common) %>%
  mutate(num_barcodes_unind = n()) %>%
  filter(num_barcodes_unind > 2) %>%
  mutate(sum_RNA_exp_UnInduced_2 = sum(normalized_fLP3_LacZ_UnInduced_RNA_2)/sum(normalized_fLP3_LacZ_UnInduced_DNA_2)) %>%
  ungroup() %>%
  select(most_common, sum_RNA_exp_UnInduced_2, num_barcodes_unind) %>%
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

# median of sum neg_control expression
med_neg_exp_ind_1 <- median(induce_exp$sum_RNA_exp_Induced_1)
med_neg_exp_ind_2 <- median(test$sum_RNA_exp_Induced_2)
med_neg_exp_un_1 <- median(test$sum_RNA_exp_UnInduced_1)
med_neg_exp_un_2 <- median(test$sum_RNA_exp_UnInduced_2)

# normalize to median or sum neg_control expression
induce_exp <- induce_exp %>%
  mutate(avg_num_ind_bc = (num_barcodes_ind.x + num_barcodes_ind.y)/2,
         avg_num_unind_bc = (num_barcodes_unind.x + num_barcodes_unind.y)/2) %>%
  mutate(normalized_RNA_exp_Induced_1 = sum_RNA_exp_Induced_1/med_neg_exp_ind_1,
         normalized_RNA_exp_Induced_2 = sum_RNA_exp_Induced_2/med_neg_exp_ind_2,
         normalized_RNA_exp_UnInduced_1 = sum_RNA_exp_UnInduced_1/med_neg_exp_un_1,
         normalized_RNA_exp_UnInduced_2 = sum_RNA_exp_UnInduced_2/med_neg_exp_un_2,
         normalized_RNA_exp_Induced_12 = (normalized_RNA_exp_Induced_1 + normalized_RNA_exp_Induced_2)/2,
         normalized_RNA_exp_UnInduced_12 = (normalized_RNA_exp_UnInduced_1 + normalized_RNA_exp_UnInduced_2)/2,
         ratio = normalized_RNA_exp_Induced_12/normalized_RNA_exp_UnInduced_12) %>%
  filter(normalized_RNA_exp_UnInduced_12 > 0) %>% # don't want variants with infinite fold-change
  select(name, variant, sum_RNA_exp_Induced_1, sum_RNA_exp_Induced_2, sum_RNA_exp_UnInduced_1, sum_RNA_exp_UnInduced_2,
         normalized_RNA_exp_Induced_1, normalized_RNA_exp_Induced_2, normalized_RNA_exp_UnInduced_1,
         normalized_RNA_exp_UnInduced_2, normalized_RNA_exp_Induced_12, normalized_RNA_exp_UnInduced_12, ratio, avg_num_ind_bc, avg_num_unind_bc)

# some housekeeping things
induce_exp$name <- gsub("shift", "", induce_exp$name)
induce_exp$name <- gsub(">LACI-", "", induce_exp$name)
induce_exp$name <- gsub("-flipped_rc", "", induce_exp$name)
induce_exp$name <- gsub("LOOP", "", induce_exp$name)
induce_exp$name <- gsub("bp", "", induce_exp$name)

backup <- induce_exp
write.table(backup, "../processed_data/LacZ_variant_exp.txt", quote = F, row.names = F)

# =========================================================================================================
# Histograms
induce_combo <- filter(induce_exp, grepl("COMBO", induce_exp$name))
induce_distal <- filter(induce_exp, grepl("DISTAL", induce_exp$name))
induce_multiple <- filter(induce_exp, grepl("MULTIPLE", induce_exp$name))
induce_steric <- filter(induce_exp, grepl("STERIC", induce_exp$name))

# COMBO LIBRARY
med_barcode_unind <- median(filter(induce_exp, grepl("COMBO", induce_exp$name))$avg_num_unind_bc)
med_barcode_ind <- median(filter(induce_exp, grepl("COMBO", induce_exp$name))$avg_num_ind_bc)
h1 <- induce_combo %>%
  gather(Condition, Count, 'avg_num_ind_bc':'avg_num_unind_bc') %>%
  ggplot(aes(Count, fill = Condition)) + 
  scale_fill_manual(values = c("#DF678C", "#3D155F"), labels = c('Induced', 'Uninduced')) +
  geom_histogram(alpha = 0.7, position = 'identity', color = 'black', binwidth = 1) +
  geom_vline(xintercept = med_barcode_unind, linetype = 'dashed') +
  annotate("text", x =15, y = 150, label = 'Median = 9', size = 4.5) +
  labs(x = 'Number of Barcodes', y = 'Count', title = 'Pcombo')

# DISTAL LIBRARY
med_barcode_unind <- median(filter(induce_exp, grepl("DISTAL", induce_exp$name))$avg_num_unind_bc)
med_barcode_ind <- median(filter(induce_exp, grepl("DISTAL", induce_exp$name))$avg_num_ind_bc)
h2 <- induce_distal %>%
  gather(Condition, Count, 'avg_num_ind_bc':'avg_num_unind_bc') %>%
  ggplot(aes(Count, fill = Condition)) + 
  scale_fill_manual(values = c("#DF678C", "#3D155F"), labels = c('Induced', 'Uninduced')) +
  geom_histogram(alpha = 0.7, position = 'identity', color = 'black', binwidth = 1) +
  geom_vline(xintercept = med_barcode_ind, linetype = 'dashed') +
  geom_vline(xintercept = med_barcode_unind, linetype = 'dashed') +
  annotate("text", x =17, y = 400, label = 'Median = 7', size = 4.5) +
  labs(x = 'Number of Barcodes', y = 'Count', title = 'Pcore')

# MULTIPLE LIBRARY
med_barcode_unind <- median(filter(induce_exp, grepl("MULTIPLE", induce_exp$name))$avg_num_unind_bc)
med_barcode_ind <- median(filter(induce_exp, grepl("MULTIPLE", induce_exp$name))$avg_num_ind_bc)
h3 <- induce_multiple %>%
  gather(Condition, Count, 'avg_num_ind_bc':'avg_num_unind_bc') %>%
  ggplot(aes(Count, fill = Condition)) + 
  scale_fill_manual(values = c("#DF678C", "#3D155F"), labels= c('Induced', 'Uninduced')) +
  geom_histogram(alpha = 0.7, position = 'identity', color = 'black', binwidth = 1) +
  geom_vline(xintercept = med_barcode_unind, linetype = 'dashed') +
  annotate("text", x =20, y = 150, label = 'Induced Median = 9', size = 4.5) +
  annotate("text", x =20, y = 130, label = 'Induced Median = 8.5', size = 4.5) +
  labs(x = 'Number of Barcodes', y = 'Count', title ='Pmultiple')

# STERIC LIBRARY
med_barcode_unind <- median(filter(induce_exp, grepl("STERIC", induce_exp$name))$avg_num_unind_bc)
med_barcode_ind <- median(filter(induce_exp, grepl("STERIC", induce_exp$name))$avg_num_ind_bc)
h4 <- induce_steric %>%
  gather(Condition, Count, 'avg_num_unind_bc':'avg_num_ind_bc') %>%
  ggplot(aes(Count, fill = Condition)) + 
  scale_fill_manual(values = c("#DF678C", "#3D155F"), labels = c('Induced', 'Uninduced')) +
  geom_histogram(alpha = 0.7, position = 'identity', color = 'black', binwidth = 1) +
  geom_vline(xintercept = med_barcode_unind, linetype = 'dashed') +
  annotate("text", x =14, y = 150, label = 'Median = 8', size = 4.5) +
  labs(x = 'Number of Barcodes', y = 'Count', title = 'Psteric')

# ======================================================================================================
# CORRELATION PLOTS

# COMBO
corr_I <- cor(induce_combo$sum_RNA_exp_Induced_1, induce_combo$sum_RNA_exp_Induced_2) #Pearson correlation
corr_UI <- cor(induce_combo$sum_RNA_exp_UnInduced_1, induce_combo$sum_RNA_exp_UnInduced_2) #Pearson correlation
# NORMALIZED INDUCED
c1a <- ggplot(induce_combo, aes(x = sum_RNA_exp_Induced_1, y = sum_RNA_exp_Induced_2)) + geom_point(alpha = .5, size = 3) +
  geom_smooth(method=lm, color = '#553C8B', se = FALSE) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_I, 3)),   parse = T, size = 5) + 
  annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Induced Biological Replicate #1') + ylab('Induced Biological Replicate #2') + labs(title = 'Pcombo') +
  geom_point(data=test, aes(sum_RNA_exp_Induced_1, sum_RNA_exp_Induced_2), color = "firebrick", size = 3, alpha = 0.5)
# NORMALIZED UNINDUCED
c1b <- ggplot(induce_combo, aes(x = sum_RNA_exp_UnInduced_1, y = sum_RNA_exp_UnInduced_2)) + geom_point(alpha = .5, size = 3) +
  geom_smooth(method=lm, color = '#553C8B', se = FALSE) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_UI, 3)),   parse = T, size = 5) + 
  annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Uninduced Biological Replicate #1') + ylab('Uninduced Biological Replicate #2') + labs(title = 'Pcombo') +
  geom_point(data=test, aes(sum_RNA_exp_UnInduced_1, sum_RNA_exp_UnInduced_2), color = "firebrick", size = 3, alpha = 0.5)

# DISTAL
corr_I <- cor(induce_distal$sum_RNA_exp_Induced_1, induce_distal$sum_RNA_exp_Induced_2) #Pearson correlation
corr_UI <- cor(induce_distal$sum_RNA_exp_UnInduced_1, induce_distal$sum_RNA_exp_UnInduced_2) #Pearson correlation
# NORMALIZED INDUCED
c2a <- ggplot(induce_distal, aes(x = sum_RNA_exp_Induced_1, y = sum_RNA_exp_Induced_2)) + geom_point(alpha = .5, size = 3) +
  geom_smooth(method=lm, color = '#553C8B', se = FALSE) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_I, 3)),   parse = T, size = 5) + 
  annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Induced Biological Replicate #1') + ylab('Induced Biological Replicate #2') + labs(title = 'Pcore') +
  geom_point(data=test, aes(sum_RNA_exp_Induced_1, sum_RNA_exp_Induced_2), color = "firebrick", size = 3, alpha = 0.5)
# NORMALIZED UNINDUCED
c2b <- ggplot(induce_distal, aes(x = sum_RNA_exp_UnInduced_1, y = sum_RNA_exp_UnInduced_2)) + geom_point(alpha = .5, size = 3) +
  geom_smooth(method=lm, color = '#553C8B', se = FALSE) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_UI, 3)),   parse = T, size = 5) + 
  annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Uninduced Biological Replicate #1') + ylab('Uninduced Biological Replicate #2') + labs(title = 'Pcore') +
  geom_point(data=test, aes(sum_RNA_exp_UnInduced_1, sum_RNA_exp_UnInduced_2), color = "firebrick", size = 3, alpha = 0.5)

# MULTIPLE
corr_I <- cor(induce_multiple$sum_RNA_exp_Induced_1, induce_multiple$sum_RNA_exp_Induced_2) #Pearson correlation
corr_UI <- cor(induce_multiple$sum_RNA_exp_UnInduced_1, induce_multiple$sum_RNA_exp_UnInduced_2) #Pearson correlation
# NORMALIZED INDUCED
c3a <- ggplot(induce_multiple, aes(x = sum_RNA_exp_Induced_1, y = sum_RNA_exp_Induced_2)) + geom_point(alpha = .5, size = 3) +
  geom_smooth(method=lm, color = '#553C8B', se = FALSE) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_I, 3)),   parse = T, size = 5) + 
  annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Induced Biological Replicate #1') + ylab('Induced Biological Replicate #2') + labs(title = 'Pmultiple') +
  geom_point(data=test, aes(sum_RNA_exp_Induced_1, sum_RNA_exp_Induced_2), color = "firebrick", size = 3, alpha = 0.5)
# NORMALIZED UNINDUCED
c3b <- ggplot(induce_multiple, aes(x = sum_RNA_exp_UnInduced_1, y = sum_RNA_exp_UnInduced_2)) + geom_point(alpha = .5, size = 3) +
  geom_smooth(method=lm, color = '#553C8B', se = FALSE) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_UI, 3)),   parse = T, size = 5) + 
  annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Uninduced Biological Replicate #1') + ylab('Uninduced Biological Replicate #2') + labs(title = 'Pmultiple') +
  geom_point(data=test, aes(sum_RNA_exp_UnInduced_1, sum_RNA_exp_UnInduced_2), color = "firebrick", size = 3, alpha = 0.5)

# STERIC
corr_I <- cor(induce_steric$sum_RNA_exp_Induced_1, induce_steric$sum_RNA_exp_Induced_2) #Pearson correlation
corr_UI <- cor(induce_steric$sum_RNA_exp_UnInduced_1, induce_steric$sum_RNA_exp_UnInduced_2) #Pearson correlation
cor.test(induce_steric$normalized_RNA_exp_UnInduced_1, induce_steric$normalized_RNA_exp_UnInduced_2)
# NORMALIZED INDUCED
c4a <- ggplot(induce_steric, aes(x = sum_RNA_exp_Induced_1, y = sum_RNA_exp_Induced_2)) + geom_point(alpha = .5, size = 3) +
  geom_smooth(method=lm, color = '#553C8B', se = FALSE) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_I, 3)),   parse = T, size = 5) + 
  annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Induced Biological Replicate #1') + ylab('Induced Biological Replicate #2') + labs(title = 'Psteric') +
  geom_point(data=test, aes(sum_RNA_exp_Induced_1, sum_RNA_exp_Induced_2), color = "firebrick", size = 3, alpha = 0.5)
# NORMALIZED UNINDUCED
c4b <- ggplot(induce_steric, aes(x = sum_RNA_exp_UnInduced_1, y = sum_RNA_exp_UnInduced_2)) + geom_point(alpha = .5, size = 3) +
  geom_smooth(method=lm, color = '#553C8B', se = FALSE) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_UI, 3)),   parse = T, size = 5) + 
  annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Uninduced Biological Replicate #1') + ylab('Uninduced Biological Replicate #2') + labs(title = 'Psteric') +
  geom_point(data=test, aes(sum_RNA_exp_UnInduced_1, sum_RNA_exp_UnInduced_2), color = "firebrick", size = 3, alpha = 0.5)

# ===========================================================================================================================================
combo <- plot_grid(h1, c1a, c1b, ncol=3)
distal <- plot_grid(h2, c2a, c2b, ncol=3)
multiple <- plot_grid(h3, c3a, c3b, ncol=3)
steric <- plot_grid(h4, c4a, c4b, ncol=3)
plot_grid(combo,distal,multiple,steric,ncol=1,labels='AUTO')
ggsave('~/Desktop/InducibleArchitectures/FreshFigs/supp_histcorr_new.pdf', height = 16, width = 17)
rm(h1, c1a, c1b, h2, c2a, c2b, h3, c3a, c3b, h4, c4a, c4b)


library(GGally)
# correlation plots
comp <- read.table('~/Desktop/newFlowComparison.txt', header = TRUE)
comp <- comp %>% select(-Name)
ggpairs(comp, columnLabels = c("Flow", "RNAseq"), lower = list(
  continuous = "smooth"))

cor(induce_exp$median_ratio, induce_exp$sum_ratio) #Pearson correlation
induce_exp %>%
  ggplot(aes(x = log10(median_ratio), y = log10(sum_ratio))) +
  geom_point(alpha = 0.75, fill = 'black')
