library(mixOmics)

## LOAD DATA

cor_all = list()

all_loads_dna = list()

all_loads_mets = list()

all_train_errors = list()


for (i in 0:4) {
  
  Xtrain_dna_name = paste0(c('XTRAIN_RNASEQ_ALL_GENES_NOREPS_SPLIT_'), i, '.csv')
  Xtrain_dna = read.csv(Xtrain_dna_name, header = TRUE, sep = ",", row.names=1)
  dim(Xtrain_dna)
  
  Xtest_dna_name = paste0(c('XTEST_RNASEQ_ALL_GENES_NOREPS_SPLIT_'), i, '.csv')
  Xtest_dna = read.csv(Xtest_dna_name, header = TRUE, sep = ",", row.names=1)
  dim(Xtest_dna)
  
  Xtrain_met_name = paste0(c('XTRAIN_METABOLOMICS_NOREPS_VT_SPLIT_'), i, '.csv')
  Xtrain_met = read.csv(Xtrain_met_name, header = TRUE, sep = ",", row.names=1)
  dim(Xtrain_met)
  
  Xtest_met_name = paste0(c('XTEST_METABOLOMICS_NOREPS_VT_SPLIT_'), i, '.csv')
  Xtest_met = read.csv(Xtest_met_name, header = TRUE, sep = ",", row.names=1)
  dim(Xtest_met)
  
  ytrain_name = paste0(c('yTRAIN_ALL_GENES_NOREPS_SPLIT_'), i, '.csv')
  ytrain = read.csv(ytrain_name, header = TRUE, sep = ",", row.names=1)
  dim(ytrain)
  
  ytest_name = paste0(c('yTEST_ALL_500_GENES_NOREPS_SPLIT_'), i, '.csv')
  ytest = read.csv(ytest_name, header = TRUE, sep = ",", row.names=1)
  dim(ytest)
  
  Y_train = ytrain[['state']]
  Y_test = ytest[['state']]
  
  train_data = list(RNASeq = Xtrain_dna,
                    metabolomics = Xtrain_met)
  
  lapply(train_data, dim)
  
  design1 = matrix(0.8, ncol = length(train_data), nrow = length(train_data), 
                   dimnames = list(names(train_data), names(train_data)))
  diag(design1) = 0 # set diagonal to 0s
  
  pls1 <- spls(train_data[["RNASeq"]], train_data[["metabolomics"]])
  
  cor1 = cor(pls1$variates$X, pls1$variates$Y)
  
  corrs = c(diag(cor1))
  
  cor_all[[i+1]] = corrs
  
  #plotVar(pls1, cutoff = 0.5, title = "RNASeq and Metabolomics", 
  #        legend = c("RNASeq", "metabolomics"), 
  #        var.names = FALSE, style = 'graphics', 
  #        pch = c(16, 17), cex = c(2,2), 
  #        col = c('coral', 'cyan3'))
  
  n_feats = list(RNASeq = c(25, 50, 100, 200, 208, 250), 
               metabolomics = c(25, 50, 100, 200, 208))
  
  tune_feats = tune.block.splsda(X = train_data, Y = Y_train, ncomp = 2, 
                                 design = design1, test.keepX = n_feats,
                                 validation = 'Mfold', folds = 10, nrepeat = 10,
                                 dist = "centroids.dist")
  
  res = tune_feats$choice.keepX # set the optimal values of features to retain
  
  keepX = res
  
  print(keepX)
  
  final.diablo.model = block.splsda(X = train_data, Y = Y_train, ncomp = 2, 
                                    design = design1, keepX = keepX)
  
  Y = final.diablo.model$Y
  
  plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median', ndisplay = 15, size.legend = 1, size.name = 1)
  plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median', ndisplay = 15, size.legend = 1, size.name = 1)
  
  loads_dna = final.diablo.model$loadings$RNASeq
  loads_mets = final.diablo.model$loadings$metabolomics
  
  all_loads_dna[[i+1]] = loads_dna
  
  all_loads_mets[[i+1]] = loads_mets
  
  data_test = list(RNASeq = Xtest_dna,
                   metabolomics = Xtest_met)
  
  perf.diablo = perf(final.diablo.model, validation = 'Mfold', 
                     M = 10, nrepeat = 10, 
                    auc = TRUE) 
  
  error = perf.diablo$WeightedVote.error.rate
  
  all_train_errors[[i+1]] = error
  
  predict.diablo = predict(final.diablo.model, newdata = data_test)
  
  
  confusion.mat1 = get.confusion_matrix(truth = Y_test,
                                       predicted = predict.diablo$WeightedVote$max.dist[,1])
  
  confusion.mat1x = get.confusion_matrix(truth = Y_test,
                                        predicted = predict.diablo$WeightedVote$max.dist[,2])
  
  
  confusion.mat2 = get.confusion_matrix(truth = Y_test,
                                        predicted = predict.diablo$WeightedVote$centroids.dist[,1])
  
  confusion.mat2x = get.confusion_matrix(truth = Y_test,
                                        predicted = predict.diablo$WeightedVote$centroids.dist[,2])
  
  confusion.mat3 = get.confusion_matrix(truth = Y_test,
                                        predicted = predict.diablo$WeightedVote$mahalanobis.dist[,1])
  
  confusion.mat3x = get.confusion_matrix(truth = Y_test,
                                        predicted = predict.diablo$WeightedVote$mahalanobis.dist[,2])
  
  print(get.BER(confusion.mat1))
  
  print(get.BER(confusion.mat1x))
  
  print(get.BER(confusion.mat2))
  
  print(get.BER(confusion.mat2x))
  
  print(get.BER(confusion.mat3))
  
  print(get.BER(confusion.mat3x))
  
  print(confusion.mat1)
  
  print(predict.diablo$WeightedVote$max.dist[,1])
  
  print(confusion.mat1x)
  
  print(predict.diablo$WeightedVote$max.dist[,2])
  
  print(confusion.mat2)
  
  print(predict.diablo$WeightedVote$centroids.dist[,1])
  
  print(confusion.mat2x)
  
  print(predict.diablo$WeightedVote$centroids.dist[,2])
  
  print(confusion.mat3)
  
  print(predict.diablo$WeightedVote$mahalanobis.dist[,1])
  
  print(confusion.mat3x)
  
  print(predict.diablo$WeightedVote$mahalanobis.dist[,2])
  
  print('-----------------------')
  
  
  
}

cor_all

# repeat for all folds

write.csv(all_loads_dna[[1]], 
          file="splits/SPLITS_DNA_LOADINGS_1.csv")

write.csv(all_loads_mets[[1]], 
          file="splits/SPLITS_METS_LOADINGS_1.csv")


all_train_errors[[1]]

# ------------------------------------------ USE DATA SPLIT 1 (2º) -------------------------------------

Xtrain_dna_name = paste0(c('XTRAIN_RNASEQ_ALL_GENES_NOREPS_SPLIT_'), 1, '.csv')
Xtrain_dna = read.csv(Xtrain_dna_name, header = TRUE, sep = ",", row.names=1)
dim(Xtrain_dna)

Xtest_dna_name = paste0(c('XTEST_RNASEQ_ALL_GENES_NOREPS_SPLIT_'), 1, '.csv')
Xtest_dna = read.csv(Xtest_dna_name, header = TRUE, sep = ",", row.names=1)
dim(Xtest_dna)

Xtrain_met_name = paste0(c('XTRAIN_METABOLOMICS_NOREPS_VT_SPLIT_'), 1, '.csv')
Xtrain_met = read.csv(Xtrain_met_name, header = TRUE, sep = ",", row.names=1)
dim(Xtrain_met)

Xtest_met_name = paste0(c('XTEST_METABOLOMICS_NOREPS_VT_SPLIT_'), 1, '.csv')
Xtest_met = read.csv(Xtest_met_name, header = TRUE, sep = ",", row.names=1)
dim(Xtest_met)

ytrain_name = paste0(c('yTRAIN_ALL_GENES_NOREPS_SPLIT_'), 1, '.csv')
ytrain = read.csv(ytrain_name, header = TRUE, sep = ",", row.names=1)
dim(ytrain)

ytest_name = paste0(c('yTEST_ALL_500_GENES_NOREPS_SPLIT_'), 1, '.csv')
ytest = read.csv(ytest_name, header = TRUE, sep = ",", row.names=1)
dim(ytest)

Y_train = ytrain[['state']]
Y_test = ytest[['state']]

train_data = list(RNASeq = Xtrain_dna,
                  metabolomics = Xtrain_met)

lapply(train_data, dim)

design1 = matrix(0.8, ncol = length(train_data), nrow = length(train_data), 
                 dimnames = list(names(train_data), names(train_data)))
diag(design1) = 0 # set diagonal to 0s

keepX = list(RNASeq = c(250,25),
             metabolomics = c(25,25))


final.diablo.model = block.splsda(X = train_data, Y = Y_train, ncomp = 2, 
                                  design = design1, keepX = keepX)

Y = final.diablo.model$Y

plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median', ndisplay = 15, size.legend = 1, size.name = 1)
plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median', ndisplay = 15, size.legend = 1, size.name = 1)

plotDiablo(final.diablo.model, ncomp = 1)
plotDiablo(final.diablo.model, ncomp = 2)

# correlations between features of the dataset  
corr = circosPlot(final.diablo.model, cutoff = 0.90, line = TRUE,
                  color.blocks= c('cyan3', 'chartreuse1'),
                  color.cor = c("chocolate3","grey20"), size.labels = 1.5, 
                  size.variables = 0.40)

# save corr between features of the dataset
write.csv(corr, 
          file="splits/RNA_METS_correlation_SPLIT1.csv")