library(mixOmics)

## LOAD DATA

Xtrain_dna = read.csv('XTRAIN_RNASEQ_MODEL_500_GENES_NOREPS.csv', header = TRUE, sep = ",", row.names=1)
dim(Xtrain_dna)

Xtest_dna = read.csv('XTEST_RNASEQ_MODEL_500_GENES_NOREPS.csv', header = TRUE, sep = ",", row.names=1)
dim(Xtest_dna)

Xtrain_flux = read.csv('XTRAIN_FLUXOMICS_500_REACTIONS.csv', header = TRUE, sep = ",", row.names=1)
dim(Xtrain_flux)

Xtest_flux = read.csv('XTEST_FLUXOMICS_500_REACTIONS.csv', header = TRUE, sep = ",", row.names=1)
dim(Xtest_flux)

ytrain = read.csv('yTRAIN_MODEL_500_GENES_NOREPS.csv', header = TRUE, sep = ",", row.names=1)
dim(ytrain)

ytest = read.csv('yTEST_MODEL_500_GENES_NOREPS.csv', header = TRUE, sep = ",", row.names=1)
dim(ytest)

Y_train = ytrain[['state']]
length(Y_train)

Y_test = ytest[['state']]
length(Y_test)

train_data = list(RNASeq = Xtrain_dna, 
                  fluxomics = Xtrain_flux)

lapply(train_data, dim) # check their dimensions


# PAIRWISE PLS (this time will select all features instead of the 25 here)


# generate pairwise PLS models
# pls1 <- spls(data[["dna"]], data[["mets"]], 
#             keepX = list.keepX, keepY = list.keepY)

pls1 <- spls(train_data[["RNASeq"]], train_data[["fluxomics"]], near.zero.var = TRUE)


plotVar(pls1, cutoff = 0.5, title = "RNASeq and Fluxomics", 
        legend = c("RNASeq", "fluxomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('coral', 'cyan3'))

cor(pls1$variates$X, pls1$variates$Y)


# design matrix: 3 options: 0.9, 0.1, and 0.5

# for square matrix filled with 0.1 ??
design1 = matrix(0.8, ncol = length(train_data), nrow = length(train_data), 
                 dimnames = list(names(train_data), names(train_data)))
diag(design1) = 0 # set diagonal to 0s

design1

design2 = matrix(0.1, ncol = length(train_data), nrow = length(train_data), 
                 dimnames = list(names(train_data), names(train_data)))
diag(design2) = 0 # set diagonal to 0s

design2

design3 = matrix(0.5, ncol = length(train_data), nrow = length(train_data), 
                 dimnames = list(names(train_data), names(train_data)))
diag(design3) = 0 # set diagonal to 0s

design3


# basic model 

# dá warning de sd is zero mas nao consigo colocar este argumento: near.zero.var = TRUE, senao perf nao corre ?

basic.diablo.model = block.splsda(X = train_data, Y = Y_train, ncomp = 5, design = design1) 

basic.diablo.model$design

# find number of components
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 5, nrepeat = 10) 

plot(perf.diablo) # plot output of tuning

perf.diablo$WeightedVote.error.rate

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
perf.diablo$choice.ncomp$WeightedVote

# find features

testn = list(RNASeq = c(25, 50, 100, 200), 
             fluxomics = c(25, 50, 100, 200))

tune_feats = tune.block.splsda(X = train_data, Y = Y_train, ncomp = 2, 
                               design = design1, test.keepX = testn,
                               validation = 'Mfold', folds = 10, nrepeat = 10,
                               dist = "centroids.dist")

res = tune_feats$choice.keepX # set the optimal values of features to retain
res

keepX = res
keepX

# final model with all features
final.diablo.model = block.splsda(X = train_data, Y = Y_train, ncomp = 2, 
                                  design = design1, keepX = keepX)

final.diablo.model$design

# see selected variables (in this case I kept it all)
selectVar(final.diablo.model, block = 'fluxomics', comp = 1)$fluxomics$name
selectVar(final.diablo.model, block = 'fluxomics', comp = 2)$fluxomics$name
selectVar(final.diablo.model, block = 'RNASeq', comp = 1)$RNASeq$name
selectVar(final.diablo.model, block = 'RNASeq', comp = 2)$RNASeq$name


final.diablo.model
Y = final.diablo.model$Y

# plotting
plotDiablo(final.diablo.model, ncomp = 1)
plotDiablo(final.diablo.model, ncomp = 2)


plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots')

# parece-me identico ao grafico anterior, vantajoso apenas quando integrar 3 omicas em vez de 2

plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'brown1'), cutoff = 0.5)

# correlations between features of the dataset  
corr = circosPlot(final.diablo.model, cutoff = 0.95, line = TRUE,
                  color.blocks= c('cyan3', 'chartreuse1'),
                  color.cor = c("chocolate3","grey20"), size.labels = 1.5, size.variables = 0.40)

network(final.diablo.model, blocks = c(1,2),
        color.node = c('cyan3', 'coral'), cutoff = 0.98)

# save corr between features of the dataset
write.csv(corr, 
          file="correlation_RNASEQ_fluxomics_featselect.csv")

# see the impact of each feature in the component 1
plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median', ndisplay = 25, size.legend = 1, size.name = 1)
plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median', ndisplay = 25, size.legend = 1, size.name = 1)

loads = final.diablo.model$loadings$RNASeq

write.csv(loads, 
          file="loading_weight.csv")

# heatmap - only after feature selection
cimDiablo(final.diablo.model, margins = c(12, 16), size.legend = 0.6)

# model evaluation
perf.diablo = perf(final.diablo.model, validation = 'Mfold', 
                   M = 10, nrepeat = 10, 
                   auc = TRUE) 

# perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate
perf.diablo$error.rate

#perf.diablo$AveragedPredict.error.rate
#perf.diablo$WeightedPredict.error.rate

perf.diablo$auc

#auroc(final.diablo.model, roc.block = c(1, 2))

# auroc por cada dataset
auc.splsda = auroc(final.diablo.model, roc.block = "RNASeq", 
                   roc.comp = 1, print = FALSE)


auroc(final.diablo.model, roc.block = "fluxomics", 
      roc.comp = 1, print = FALSE)



# check performance TEST DATA

data_test = list(RNASeq = Xtest_dna,
                 fluxomics = Xtest_flux)

predict.diablo = predict(final.diablo.model, newdata = data_test)

predict.diablo$MajorityVote
predict.diablo$WeightedVote

confusion.mat = get.confusion_matrix(truth = Y_test,
                                     predicted = predict.diablo$WeightedVote$centroids.dist[,1])
confusion.mat

get.BER(confusion.mat)

Y_test
