library(dplyr)
library(tidyr)

induce_exp <- read.table("~/Desktop/InducibleArchitectures/processed_data/induce_exp(updated).txt", header = T)
induce_combo <- subset(induce_exp, grepl("COMBO", induce_exp$name)) %>% separate(col = 'name', into = c("Library", "Distal", "Min_35", "Min_10", "Proximal"), sep = "-", remove = F)

modeling_combo <- induce_combo %>%
  select(Distal, Min_35, Min_10, Proximal, normalized_RNA_exp_UnInduced_12, normalized_RNA_exp_Induced_12, ratio)

#modeling_combo$Distal <- paste0('"', modeling_combo$Distal, '",')
#modeling_combo$Min_35 <- paste0('"', modeling_combo$Min_35, '",')
#modeling_combo$Min_10 <- paste0('"', modeling_combo$Min_10, '",')
#modeling_combo$Proximal <- paste0('"', modeling_combo$Proximal, '",')
#modeling_combo$normalized_RNA_exp_UnInduced_12 <- paste0(modeling_combo$normalized_RNA_exp_UnInduced_12, "`,")
#modeling_combo$normalized_RNA_exp_Induced_12 <- paste0(modeling_combo$normalized_RNA_exp_Induced_12, "`")

#modeling_combo$normalized_RNA_exp_UnInduced_12 <- paste0("'", modeling_combo$normalized_RNA_exp_UnInduced_12, "',")
#modeling_combo$normalized_RNA_exp_Induced_12 <- paste0("'", modeling_combo$normalized_RNA_exp_Induced_12, "'")

write.table(modeling_combo, "~/Desktop/InducibleArchitectures/processed_data/induce_combo.txt", row.names = FALSE, col.names = FALSE, sep = '\t')

nulldistal_combo <- induce_combo %>%
  filter(Distal == 'lacOscram') %>%
  select(Distal, Min_35, Min_10, Proximal, normalized_RNA_exp_UnInduced_12, normalized_RNA_exp_Induced_12) 
write.csv(nulldistal_combo, "~/Desktop/InducibleArchitectures/processed_data/induce_combo_nulldistal.csv", row.names = FALSE)

