barcode_stats_loop <- read.table("~/Desktop/InducibleArchitectures/raw_data/LacZ_mapping_barcode_statistics.txt", header = T)
mapped_barcodes <- barcode_stats_loop[!is.na(barcode_stats_loop$most_common),] # removed unmapped
Compare_barcode_Reps <- left_join(mapped_barcodes, fLP3_loop, by ='barcode')

Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0 #Make all NA values = 0
temp <- Compare_barcode_Reps # save it in temp dataframe 

# FUNCTION FOR CREATING DATAFRAME MADE WITH SUMMED EXP, num_barcodes > 2 filter
sum_var <- function(df) {
  uninduced <- df %>%
    select(most_common, normalized_fLP3_LacZ_UnInduced_DNA_1, normalized_fLP3_LacZ_UnInduced_DNA_2, normalized_fLP3_LacZ_UnInduced_RNA_1,
           normalized_fLP3_LacZ_UnInduced_RNA_2) %>%
    filter(normalized_fLP3_LacZ_UnInduced_DNA_1 > 0 & 
             normalized_fLP3_LacZ_UnInduced_DNA_2 > 0) %>%
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
    filter(normalized_fLP3_LacZ_Induced_DNA_1 > 0 & 
             normalized_fLP3_LacZ_Induced_DNA_2 > 0) %>%
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

table <- temp %>% filter(most_common == 'GCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGAATGTAAGTTAGGAAAATTTTTTTTCAAAAGTAGCAATTGTGAGCGGATAACAATTTGGTATAATGTGTGGACACACAGGAAACAGCTATAATTGTGAGCGGATAACAATTGACC') %>%
  #filter(normalized_fLP3_LacZ_UnInduced_DNA_1 > 0 & normalized_fLP3_LacZ_UnInduced_DNA_2 > 0 &
  #       normalized_fLP3_LacZ_Induced_DNA_1 > 0 & normalized_fLP3_LacZ_Induced_DNA_2 > 0) %>%
  #mutate(DNA_UnInduced_avg = mean(sum(normalized_fLP3_LacZ_UnInduced_DNA_1), sum(normalized_fLP3_LacZ_UnInduced_DNA_2))) %>%
  mutate(RNA_exp_UnInduced_1 = normalized_fLP3_LacZ_UnInduced_RNA_1/normalized_fLP3_LacZ_UnInduced_DNA_1,
         RNA_exp_UnInduced_2 = normalized_fLP3_LacZ_UnInduced_RNA_2/normalized_fLP3_LacZ_UnInduced_DNA_2) %>%
  #mutate(DNA_Induced_avg = mean(sum(normalized_fLP3_LacZ_Induced_DNA_1), sum(normalized_fLP3_LacZ_Induced_DNA_2))) %>%
  mutate(RNA_exp_Induced_1 = normalized_fLP3_LacZ_Induced_RNA_1/normalized_fLP3_LacZ_Induced_DNA_1,
         RNA_exp_Induced_2 = normalized_fLP3_LacZ_Induced_RNA_2/normalized_fLP3_LacZ_Induced_DNA_2)
  #select(most_common, RNA_exp_UnInduced_1, RNA_exp_UnInduced_2, RNA_exp_Induced_1, RNA_exp_Induced_2)

test <- subset(induce_exp, grepl("neg_control", induce_exp$name)) 

med_neg_exp_ind_1 <- median(test$RNA_exp_Induced_1)
med_neg_exp_ind_2 <- median(test$RNA_exp_Induced_2)
med_neg_exp_un_1 <- median(test$RNA_exp_UnInduced_1)
med_neg_exp_un_2 <- median(test$RNA_exp_UnInduced_2)

table <- table %>%
  mutate(normalized_RNA_exp_Induced_1 = RNA_exp_Induced_1/med_neg_exp_ind_1,
         normalized_RNA_exp_Induced_2 = RNA_exp_Induced_2/med_neg_exp_ind_2,
         normalized_RNA_exp_UnInduced_1 = RNA_exp_UnInduced_1/med_neg_exp_un_1,
         normalized_RNA_exp_UnInduced_2 = RNA_exp_UnInduced_2/med_neg_exp_un_2,
         normalized_RNA_exp_Induced_12 = (normalized_RNA_exp_Induced_1 + normalized_RNA_exp_Induced_2)/2,
         normalized_RNA_exp_UnInduced_12 = (normalized_RNA_exp_UnInduced_1 + normalized_RNA_exp_UnInduced_2)/2,
         ratio = normalized_RNA_exp_Induced_12/normalized_RNA_exp_UnInduced_12) %>%
  select(most_common, normalized_RNA_exp_Induced_1, normalized_RNA_exp_Induced_2, normalized_RNA_exp_UnInduced_1,
         normalized_RNA_exp_UnInduced_2, normalized_RNA_exp_Induced_12, normalized_RNA_exp_UnInduced_12, ratio)
write.table(table, '~/Desktop/S1_table.txt', sep = '\t', row.names = FALSE)  
  
# ==============================
# C1
table <- temp %>% filter(most_common == 'ATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGAAATTGTGAGCGCTCACAATTTTCATTAGGCACCCCAGGCTTTGACATTGTGAGCGGATAACAATATAATGTGTGGACACACAGGAAACAGCTATGACCATGATT') %>%
  filter(normalized_fLP3_LacZ_UnInduced_DNA_1 > 0 & normalized_fLP3_LacZ_UnInduced_DNA_2 > 0) %>%
  filter(normalized_fLP3_LacZ_Induced_DNA_1 > 0 & normalized_fLP3_LacZ_Induced_DNA_2 > 0) %>%
  mutate(DNA_UnInduced_avg = mean(sum(normalized_fLP3_LacZ_UnInduced_DNA_1), sum(normalized_fLP3_LacZ_UnInduced_DNA_2))) %>%
  mutate(RNA_exp_UnInduced_1 = normalized_fLP3_LacZ_UnInduced_RNA_1/normalized_fLP3_LacZ_UnInduced_DNA_1,
         RNA_exp_UnInduced_2 = normalized_fLP3_LacZ_UnInduced_RNA_2/normalized_fLP3_LacZ_UnInduced_DNA_2) %>%
  #mutate(DNA_Induced_avg = mean(sum(normalized_fLP3_LacZ_Induced_DNA_1), sum(normalized_fLP3_LacZ_Induced_DNA_2))) %>%
  mutate(RNA_exp_Induced_1 = normalized_fLP3_LacZ_Induced_RNA_1/normalized_fLP3_LacZ_Induced_DNA_1,
         RNA_exp_Induced_2 = normalized_fLP3_LacZ_Induced_RNA_2/normalized_fLP3_LacZ_Induced_DNA_2) %>%
  select(most_common, RNA_exp_UnInduced_1, RNA_exp_UnInduced_2, RNA_exp_Induced_1, RNA_exp_Induced_2)

test <- subset(induce_exp, grepl("neg_control", induce_exp$name)) 

med_neg_exp_ind_1 <- median(test$RNA_exp_Induced_1)
med_neg_exp_ind_2 <- median(test$RNA_exp_Induced_2)
med_neg_exp_un_1 <- median(test$RNA_exp_UnInduced_1)
med_neg_exp_un_2 <- median(test$RNA_exp_UnInduced_2)

table <- table %>%
  mutate(normalized_RNA_exp_Induced_1 = RNA_exp_Induced_1/med_neg_exp_ind_1,
         normalized_RNA_exp_Induced_2 = RNA_exp_Induced_2/med_neg_exp_ind_2,
         normalized_RNA_exp_UnInduced_1 = RNA_exp_UnInduced_1/med_neg_exp_un_1,
         normalized_RNA_exp_UnInduced_2 = RNA_exp_UnInduced_2/med_neg_exp_un_2,
         normalized_RNA_exp_Induced_12 = (normalized_RNA_exp_Induced_1 + normalized_RNA_exp_Induced_2)/2,
         normalized_RNA_exp_UnInduced_12 = (normalized_RNA_exp_UnInduced_1 + normalized_RNA_exp_UnInduced_2)/2,
         ratio = normalized_RNA_exp_Induced_12/normalized_RNA_exp_UnInduced_12) %>%
  select(most_common, normalized_RNA_exp_Induced_1, normalized_RNA_exp_Induced_2, normalized_RNA_exp_UnInduced_1,
         normalized_RNA_exp_UnInduced_2, normalized_RNA_exp_Induced_12, normalized_RNA_exp_UnInduced_12, ratio)

