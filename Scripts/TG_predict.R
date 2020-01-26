library(dplyr)
library(glmnet)
library(ggplot2)
library(viridis)
require(cowplot)
library(tidyr)
library(bnstruct)
library(data.table)
library(ggsignif)
library(Hmisc)

# read data in ======================================
hmdp_protein = read.table('~/Desktop/HMDP/raw_data/new_data/proteomics.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
hmdp_lipid = read.table('~/Desktop/HMDP/raw_data/new_data/HMDP_liver_lipidomics_allstrains_individual_lipid_species_for_TV_plus_identifiers.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
protein_names = read.csv('~/Desktop/HMDP/raw_data/new_data/protein_mapping.csv', header = FALSE, stringsAsFactors = FALSE)
lipid_names = read.table('~/Desktop/HMDP/raw_data/new_data/lipid_colNames.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
# formatting data ======================================
hmdp_protein[hmdp_protein=="NaN"] <- NA
protein_input = hmdp_protein[ ,colSums(is.na(hmdp_protein)) <= 145] # remove any column with NA data, ALL replicates must have measurements
protein_input <- cbind(protein_input[,c(1,2)], as.data.frame(knn.impute(as.matrix(as.data.frame(lapply(protein_input[,-c(1,2)], as.numeric))))))

lipid_input = as.data.frame(t(hmdp_lipid[,-1]))
lipid_input = setnames(lipid_input, old = colnames(lipid_input)[1:313], new = as.character(lipid_names$V1)) 

totalTG = lipid_input %>% select(1, starts_with("TG ")) %>%
  gather(Lipid, Value, 2:43) %>%
  group_by(Mouse_ID) %>%
  mutate(totalTG = sum(as.numeric(Value))) %>%
  ungroup() %>% 
  select(-Lipid, -Value) %>%
  distinct()

lasso_input = left_join(totalTG, protein_input, by = 'Mouse_ID') %>%
  select(-Mouse_ID, -Strain) %>% data.frame(stringsAsFactors = FALSE)

# exploratory analyses, heatmap of TG features, see which ones might be interesting to study ===============
lipid_input %>% select(1, starts_with("TG.")) %>%
  gather(Lipid, Value, 2:43) %>%
  group_by(Mouse_ID) %>%
  mutate(totalTG = sum(as.numeric(Value))) %>%
  ungroup() %>% 
  select(-Lipid, -Value) %>%
  distinct()

TG_features = lipid_input %>% select(starts_with("TG."))

rcorr(as.matrix(TG_features), type = 'spearman')$r %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = 'Feature1') %>%
  gather(Feature2, Corr, 2:43) %>%
  ggplot(aes(x = Feature1, y = Feature2)) + geom_tile(aes(fill = Corr), alpha = 0.9, color = 'black', size = 0.1) + 
  scale_fill_distiller(palette = 'RdYlBu', limits = c(-1,1)) + labs(title = 'Correlations between TG species') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10))
ggsave('~/Desktop/HMDP/freshFigs/TG_correlation.pdf', height = 8, width = 10)

# run regression ========================================
runLassoRegression <- function(input_data, numFolds, numRep, numPerm, lambda_idx) {
  # things to initialize =================================================================
  classification_table <- data.frame("condition" = character(),
                                     "accuracy" = double(), 
                                     stringsAsFactors=FALSE)
  # get p value of overall model with this function
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  # =======================================================================================
  # Regression part
  for(j in 1:numRep) {
    folds <- cut(seq(1,nrow(input_data)), breaks=numFolds, labels=FALSE) %>% sample(., length(.), replace= F)
    for(i in 1:numFolds) {
      foldNum = i # to keep track of where we are
      # Initialize folds, training and test datasets
      testIdx <- which(folds == i, arr.ind=TRUE)
      trainIdx <- which(folds != i, arr.ind= TRUE)
      testData <-  input_data[testIdx,]
      trainData <- input_data[-testIdx,]
      features <- as.matrix(trainData[,-1])
      class <- as.vector(trainData[,1])
      # Perform cross-validation and extract betas from LASSO
      numFeatures = 0
        cvfit = cv.glmnet(features, class, alpha = 1, intercept = TRUE, grouped = FALSE)
        betas = as.data.frame(as.matrix(coef(cvfit, s = cvfit$lambda[lambda_idx]))) %>%
          tibble::rownames_to_column(var = 'feature')
        names(betas) <- c('feature', 'beta')
        betas <- betas %>% dplyr::filter(beta != 0, feature != '(Intercept)')
        feature_list <- as.vector(betas$feature)
        print(feature_list)
        numFeatures <- length(feature_list)
        print(paste('Number features selected: ', numFeatures))
      
      # linear model significance
      if(numFeatures == 1) {
        idx <- match(feature_list, colnames(trainData))
        subsetTrainData <- as.data.frame(trainData[,idx])
        subsetTrainData <- cbind(trainData[,1], subsetTrainData)
        colnames(subsetTrainData)[1] <- "Response"
        colnames(subsetTrainData)[2] <- feature_list
      } else {
        idx <- match(feature_list, colnames(trainData))
        subsetTrainData <- trainData[,idx]
        subsetTrainData <- cbind(trainData[,1], subsetTrainData)
        colnames(subsetTrainData)[1] <- "Response"
      }

      lmfit = lm(Response~., data = subsetTrainData)
      testData$predicted_TG <- predict(lmfit, newdata=testData) #predict total triglyceride of train data
      actual_predict = cor.test(testData$predicted_TG, testData$totalTG, method = 'pearson')$p.value
      
      # significance of model trained on actual data
      classification_table[nrow(classification_table) + 1, ] <- c("actual", -log10(actual_predict))
      
      # PERMUTED AND RANDOM DATA =====================================================================
      print(paste('Round', foldNum,'/5 in permutation testing for repetition', j,'.'))
      if(numPerm > 0) {
        for(repetition in 1:numPerm) {
          perm_class = sample(class)
          perm_numFeatures = 0
          # internal check to make sure the response vector is not constant 
          perm_cvfit = cv.glmnet(features, perm_class, alpha = 1, intercept = TRUE)
          all_betas = as.data.frame(as.matrix(coef(perm_cvfit, s = perm_cvfit$lambda[lambda_idx]))) %>% 
            tibble::rownames_to_column(var = 'feature')
          names(all_betas) <- c('feature', 'beta')
          subset_betas <- all_betas %>% filter(beta != 0, feature != '(Intercept)')
          perm_feature_list <- as.vector(subset_betas$feature)
          perm_numFeatures <- length(perm_feature_list)
          if(perm_numFeatures < 1) {
            rand_betas <- all_betas %>% filter(feature != '(Intercept)') %>%
              sample_n(2, replace = FALSE)
            perm_feature_list <- as.vector(rand_betas$feature)
            perm_numFeatures <- length(perm_feature_list)
          }

          # PERMUTED linear model significance
          if(perm_numFeatures == 1) {
            perm_idx <- match(perm_feature_list, colnames(trainData))
            perm_subsetTrainData <- as.data.frame(trainData[,perm_idx])
            perm_subsetTrainData <- cbind(trainData[,1], perm_subsetTrainData)
            colnames(perm_subsetTrainData)[1] <- "Response"
            colnames(perm_subsetTrainData)[2] <- perm_feature_list
          } else {
            perm_idx <- match(perm_feature_list, colnames(trainData))
            perm_subsetTrainData <- trainData[,perm_idx] 
            perm_subsetTrainData <- cbind(trainData[,1], perm_subsetTrainData)
            colnames(perm_subsetTrainData)[1] <- "Response"
          }

          perm_lmfit = lm(Response~., data = perm_subsetTrainData)
          testData$perm_predicted_TG <- predict(perm_lmfit, newdata=testData) #predict total triglyceride of train data
          perm_predict = cor.test(testData$perm_predicted_TG, testData$totalTG, method = 'pearson')$p.value
          
          # significance of model trained on permuted data
          classification_table[nrow(classification_table) + 1, ] <- c("permuted", -log10(perm_predict))
          
          # Random part
          random_feat <- as.data.frame(colnames(features)) %>%
            sample_n(numFeatures, replace = FALSE)
          names(random_feat)[1] = 'feature'
          rand_feature_list <- as.vector(random_feat$feature)
          
          # RANDOM linear model significance
          if(numFeatures == 1) {
            rand_idx <- match(rand_feature_list, colnames(trainData))
            rand_subsetTrainData <- as.data.frame(trainData[,rand_idx])
            rand_subsetTrainData <- cbind(trainData[,1], rand_subsetTrainData)
            colnames(rand_subsetTrainData)[1] <- "Response"
            colnames(rand_subsetTrainData)[2] <-rand_feature_list
          } else {
            rand_idx <- match(rand_feature_list, colnames(trainData))
            rand_subsetTrainData <- trainData[,rand_idx] 
            rand_subsetTrainData <- cbind(trainData[,1], rand_subsetTrainData)
            colnames(rand_subsetTrainData)[1] <- "Response"
          }
          
          rand_lmfit = lm(Response~., data = rand_subsetTrainData)
          testData$rand_predicted_TG <- predict(rand_lmfit, newdata=testData) #predict total triglyceride of train data
          rand_predict = cor.test(testData$rand_predicted_TG, testData$totalTG, method = 'pearson')$p.value
          
          # significance of model trained on permuted data
          classification_table[nrow(classification_table) + 1, ] <- c("random", -log10(rand_predict))
        }
      }
    }
  }
  return(classification_table)
}
set.seed(123)
class_table <- runLassoRegression(lasso_input, 5, 10, 1, 10)
med_actual <- median(as.numeric(filter(class_table, condition == 'actual')$accuracy))
med_permute <- median(as.numeric(filter(class_table, condition == 'permuted')$accuracy))
med_random <- median(as.numeric(filter(class_table, condition == 'random')$accuracy))

class_table %>%
  ggplot(aes(x = condition, y = as.numeric(accuracy))) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = condition), pch = 21, color = 'black', size = 3.5, alpha = 1) +
  scale_fill_manual(values = c("#800075", "#6CA6CD", "#ffdf00")) +
  #scale_fill_viridis(discrete = TRUE) +
  labs(x = 'Condition', y = 'Predictability') +
  geom_signif(comparisons = list(c("actual", "permuted"), c("actual", "random")), 
              y_position= c(25.5, 23.5), map_signif_level=function(p)sprintf("p = %.2g", p), test="wilcox.test") + 
  geom_hline(yintercept = med_actual, linetype = 'dashed', color = 'red', size = 0.75) +
  geom_hline(yintercept = med_permute, linetype = 'dashed', color = 'black', size = 0.75) +
  geom_hline(yintercept = med_random, linetype = 'dashed', color = 'black', size = 0.75) +
  theme(legend.position = 'none') +
  ggtitle('Total TG Prediction')
ggsave('~/Desktop/HMDP/freshFigs/TG_model.png', height = 5, width = 5)

runFinalLassoModel <- function(input_data, lambda_idx) {
  set.seed(123)
  final_cv = cv.glmnet(as.matrix(input_data[,-1]), as.vector(input_data[,1]), nfolds = 5, alpha = 1, intercept = TRUE)
  final_betas = as.data.frame(as.matrix(coef(final_cv, s = final_cv$lambda[lambda_idx]))) %>% 
    tibble::rownames_to_column(var = 'feature')
  names(final_betas) <- c('feature', 'beta')
  final_betas <- final_betas %>% filter(beta != 0, feature != '(Intercept)')
  final_feature_list <- as.vector(final_betas$feature)
  return(final_feature_list)
}
selectFeat <- runFinalLassoModel(lasso_input, 10)
select_proteins = protein_names %>% filter(V1 %in% selectFeat)

lasso_input %>%
  ggplot(aes(x = totalTG, y = `protein182`)) + geom_point(alpha = 0.5, size = 3) + 
  geom_smooth(method=lm, color = '#553C8B', se = FALSE)
ggsave('~/Desktop/HMDP/freshFigs/PLIN2_TG.pdf', height = 5, width = 6)
# correlation network analysis =====================================
library(Hmisc)
options(digits=10)
# Function to remove duplicate values
killDuplicates <- function(x) {
  cols = c(1,2)
  temp <- x[,cols]
  for (i in 1:nrow(x)){
    temp[i, ] = sort(x[i,cols])
  }
  x <- as.data.frame(x[!duplicated(temp),])
  return(x)
}

corr <- rcorr(as.matrix(lasso_input[,-1]), type = 'spearman')$r %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = 'Feature1') %>%
  gather(Feature2, Corr, 2:2920) %>%
  filter(Feature1 != Feature2) # filter all same feature correlations

# calculate P-value
pval <-rcorr(as.matrix(lasso_input[,-1]), type = 'spearman')$P %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = 'Feature1') %>%
  gather(Feature2, Pval, 2:2920) %>%
  mutate(FDR = p.adjust(Pval, method = 'BH')) %>%
  dplyr::select(Feature1, Feature2, FDR) %>%
  filter(Feature1 != Feature2) # filter all same feature correlations

# we want the edges of our primary and secondary features (gives us all primary-primary, primary-secondary edges)
#options(digits=20)
full_corr <- left_join(corr, pval, by = c('Feature1','Feature2')) %>%
  filter(Corr < -0.4 | Corr > 0.4) %>%
  filter(FDR < 0.05) %>%
  filter(Feature1 %in% selectFeat | Feature2 %in% selectFeat) %>%
  killDuplicates()

all_feats <- unique(append(full_corr$Feature1, full_corr$Feature2)) # gives us all primary and secondary featurs

# Will select all edges between primary-secondary, primary-primary, and 
# secondary-secondary interactions that pass threshold.
cyto_corr <- left_join(corr, pval, by = c('Feature1','Feature2')) %>%
  filter(Corr < -0.4 | Corr > 0.4) %>%
  filter(FDR < 0.05) %>% 
  filter((Feature1 %in% all_feats & Feature2 %in% all_feats)) %>%
  killDuplicates()

# sanity check: should equal number nodes in full_corr
nodes <- unique(c(as.character(cyto_corr$Feature1), as.character(cyto_corr$Feature2)))
network_proteins = protein_names %>% filter(V1 %in% nodes)

# network to cytoscape
write.csv(cyto_corr, '~/Desktop/HMDP/processed_data/CE_cytoscape_net.csv')

