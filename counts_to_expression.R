setwd("~/Desktop/induce_library/raw_data/")

#install.packages("ggjoy")
#install.packages("RInside") #THis is for ggjoy
#install.packages("purrr") #This is for ggjoy
#install.packages("plotly")

#theme_update(text=(theme_grey()+theme(text = element_text(family = "sans")))$text )
library(plyr)
library("ggplot2")
library("dplyr")
library("wesanderson")
names(wes_palettes)
library("tidyr")
library("ggjoy")
library("reshape2")
require(cowplot)
library('RInside')
library("stringr")
library("plotly")


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

names(fLP3_loop) = sub(".txt","", names(fLP3_loop)) #rename all colummns that were named after text ile temp <- Loop_Data[c(-1,-5)]

rm(list = c(filelist))
rm(x)

# Read in all barcodes from barcode stats file
barcode_stats_loop <- read.table("LacZ_mapping_barcode_statistics.txt", header = T)
mapped_barcodes <- barcode_stats_loop[!is.na(barcode_stats_loop$most_common),] # removed unmapped
Compare_barcode_Reps <- left_join(mapped_barcodes, fLP3_loop, by ='barcode')

Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0 #Make all NA values = 0
temp <- Compare_barcode_Reps # save it in temp dataframe 

# FUNCTION FOR CREATING DATAFRAME MADE WITH SUMMED EXP, num_barcodes > 2

sum_var <- function(df) {
  uninduced <- df %>%
    select(most_common, normalized_fLP3_LacZ_UnInduced_DNA_1, normalized_fLP3_LacZ_UnInduced_DNA_2, normalized_fLP3_LacZ_UnInduced_RNA_1,
           normalized_fLP3_LacZ_UnInduced_RNA_2) %>%
    filter(normalized_fLP3_LacZ_UnInduced_DNA_1 > 0 & 
             normalized_fLP3_LacZ_UnInduced_DNA_2 > 0) %>%
    group_by(most_common) %>%
    mutate(num_barcodes = n()) %>%
    filter(num_barcodes > 2) %>%
    mutate(DNA_UnInduced_avg = mean(sum(normalized_fLP3_LacZ_UnInduced_DNA_1), sum(normalized_fLP3_LacZ_UnInduced_DNA_2))) %>%
    mutate(RNA_exp_UnInduced_1 = sum(normalized_fLP3_LacZ_UnInduced_RNA_1)/(DNA_UnInduced_avg),
           RNA_exp_UnInduced_2 = sum(normalized_fLP3_LacZ_UnInduced_RNA_2)/(DNA_UnInduced_avg)) %>%
    ungroup() %>%
    select(-normalized_fLP3_LacZ_UnInduced_DNA_1, -normalized_fLP3_LacZ_UnInduced_DNA_2, -normalized_fLP3_LacZ_UnInduced_RNA_1, -normalized_fLP3_LacZ_UnInduced_RNA_2) %>%
    distinct()
  
  induced <- df %>%
    select(most_common, normalized_fLP3_LacZ_Induced_DNA_1, normalized_fLP3_LacZ_Induced_DNA_2, normalized_fLP3_LacZ_Induced_RNA_1,
           normalized_fLP3_LacZ_Induced_RNA_2) %>%
    filter(normalized_fLP3_LacZ_Induced_DNA_1 > 0 & 
             normalized_fLP3_LacZ_Induced_DNA_2 > 0) %>%
    group_by(most_common) %>%
    mutate(num_barcodes = n()) %>%
    filter(num_barcodes > 2) %>%
    mutate(DNA_Induced_avg = mean(sum(normalized_fLP3_LacZ_Induced_DNA_1), sum(normalized_fLP3_LacZ_Induced_DNA_2))) %>%
    mutate(RNA_exp_Induced_1 = sum(normalized_fLP3_LacZ_Induced_RNA_1)/(DNA_Induced_avg),
           RNA_exp_Induced_2 = sum(normalized_fLP3_LacZ_Induced_RNA_2)/(DNA_Induced_avg)) %>%
    ungroup() %>%
    select(-normalized_fLP3_LacZ_Induced_DNA_1, -normalized_fLP3_LacZ_Induced_DNA_2, -normalized_fLP3_LacZ_Induced_RNA_1, -normalized_fLP3_LacZ_Induced_RNA_2) %>%
    distinct()
  
  inner_join(induced, uninduced, by = 'most_common', suffix = c('_ind', '_unind'))
}

mean_exp_master <- sum_var(temp) 

# read in variant file to get all names of promoters 
variant_stats <- read.table("LacZ_mapping_variant_statistics.txt", header = T, fill = T)

# create dataframe with all variant names and filter out expression > 0
induce_exp <- left_join(variant_stats, mean_exp_master, by = c("variant" = "most_common")) %>% na.omit(.) %>% filter(RNA_exp_Induced_1 > 0 & RNA_exp_Induced_2 > 0 & RNA_exp_UnInduced_1 > 0 & RNA_exp_UnInduced_2 > 0) %>% select("variant", "name", "RNA_exp_Induced_1", "RNA_exp_Induced_2", "RNA_exp_UnInduced_1", "RNA_exp_UnInduced_2")

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

#expression for each sample is normalized first by dividing by the median of the negative control expression in that corresponding sample. These normalized expressions are then averaged across the two replicates for Induced and UnInduced. Ratio is determined by dividing normalized Induced expression by normalized UnInduced expression. 
induce_exp <- induce_exp %>%
  mutate(normalized_RNA_exp_Induced_1 = RNA_exp_Induced_1/med_neg_exp_ind_1,
         normalized_RNA_exp_Induced_2 = RNA_exp_Induced_2/med_neg_exp_ind_2,
         normalized_RNA_exp_UnInduced_1 = RNA_exp_UnInduced_1/med_neg_exp_un_1,
         normalized_RNA_exp_UnInduced_2 = RNA_exp_UnInduced_2/med_neg_exp_un_2,
         normalized_RNA_exp_Induced_12 = (normalized_RNA_exp_Induced_1 + normalized_RNA_exp_Induced_2)/2,
         normalized_RNA_exp_UnInduced_12 = (normalized_RNA_exp_UnInduced_1 + normalized_RNA_exp_UnInduced_2)/2,
         ratio = normalized_RNA_exp_Induced_12/normalized_RNA_exp_UnInduced_12)

# correlation between Induced technical replicates
corr_I <- cor(induce_exp$RNA_exp_Induced_1, induce_exp$RNA_exp_Induced_2) #Pearson correlation
# correlation between UnInduced technical replicates
corr_UI <- cor(induce_exp$RNA_exp_UnInduced_1, induce_exp$RNA_exp_UnInduced_2) #Pearson correlation

# NORMALIZED INDUCED

ggplot(induce_exp, aes(normalized_RNA_exp_UnInduced_1, normalized_RNA_exp_UnInduced_2)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_I, 3)),   parse = T, size = 5) + annotation_logticks() + scale_x_log10(limits = c(0.01,300), breaks = c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + xlab('Normalized Induced 1') + ylab('Normalized Induced 2') + ggtitle('Comparing Normalized Promoter Expression Between Technical Replicates (Induced)') + geom_point(data=test, aes(normalized_RNA_exp_UnInduced_1, normalized_RNA_exp_UnInduced_2), color = "red")
ggsave('InducePlot.pdf', width = 10, height = 6)

# NORMALIZED UNINDUCED

ggplot(induce_exp, aes(normalized_RNA_exp_UnInduced_1, normalized_RNA_exp_UnInduced_2)) + geom_point(alpha = .2) +
  geom_smooth(method=lm) + annotate("text", x =.5, y = 30, label = paste('r==', signif(corr_UI, 3)),   parse = T,  size = 5) + annotation_logticks() + scale_x_log10(breaks = c(1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + xlab('UnInduced 1') + ylab('UnInduced 2') + ggtitle('Comparing Normalized Promoter Expression Between Technical Replicates (UnInduced)') + geom_point(data=test, aes(normalized_RNA_exp_UnInduced_1, normalized_RNA_exp_UnInduced_2), color = "red") 
ggsave('UnInducePlot.pdf', width = 10, height = 6)


induce_exp$name <- gsub("shift", "", induce_exp$name)
induce_exp$name <- gsub(">LACI-", "", induce_exp$name)
induce_exp$name <- gsub("-flipped_rc", "", induce_exp$name)
induce_exp$name <- gsub("LOOP", "", induce_exp$name)
induce_exp$name <- gsub("bp", "", induce_exp$name)

backup <- induce_exp
write.table(backup, "../processed_data/induce_exp.txt", quote = F, row.names = F)
#induce_exp <- read.table("./induce_exp.txt", header = T)


# check induced expression for Osym + Osym (maybe it likes being there so much it can't get expressed)
# change everything to log2




