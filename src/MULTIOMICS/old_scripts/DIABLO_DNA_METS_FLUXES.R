library(mixOmics)

## LOAD DATA

Xtrain_dna = read.csv('XTRAIN_RNASEQ_ALL_500_GENES_NOREPS.csv', header = TRUE, sep = ",", row.names=1)
dim(Xtrain_dna)

Xtest_dna = read.csv('XTEST_RNASEQ_ALL_500_GENES_NOREPS.csv', header = TRUE, sep = ",", row.names=1)
dim(Xtest_dna)

Xtrain_met = read.csv('XTRAIN_METABOLOMICS_NOREPS_VT.csv', header = TRUE, sep = ",", row.names=1)
dim(Xtrain_met)

Xtest_met = read.csv('XTEST_METABOLOMICS_NOREPS_VT.csv', header = TRUE, sep = ",", row.names=1)
dim(Xtest_met)

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
                  metabolomics = Xtrain_met,
                  fluxomics = Xtrain_flux)

lapply(train_data, dim) # check their dimensions


# PAIRWISE PLS (apenas entre fluxomics e mets porque é o que falta)

pls1 <- spls(train_data[["metabolomics"]], train_data[["fluxomics"]])

plotVar(pls1, cutoff = 0.5, title = "Metabolomics and Fluxomics", 
        legend = c("Metabolomics", "Fluxomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('coral', 'cyan3'))

cor(pls1$variates$X, pls1$variates$Y)

# DESIGN MATRIX: all datasets need to have the same value or can I adjust? what makes sense?

design1 = matrix(0.5, ncol = length(train_data), nrow = length(train_data), 
                 dimnames = list(names(train_data), names(train_data)))
diag(design1) = 0 # set diagonal to 0s

design1

# basic model 

basic.diablo.model = block.splsda(X = train_data, Y = Y_train, ncomp = 5, design = design1) 


# find number of components
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 5, nrepeat = 10) 

plot(perf.diablo) # plot output of tuning

perf.diablo$WeightedVote.error.rate

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
perf.diablo$choice.ncomp$WeightedVote


# find features

testn = list(RNASeq = c(25, 50, 100),
             metabolomics = c(25, 50, 100),
             fluxomics = c(25, 50, 100))

tune_feats = tune.block.splsda(X = train_data, Y = Y_train, ncomp = 2, 
                               design = design1, test.keepX = testn,
                               validation = 'Mfold', folds = 10, nrepeat = 10,
                               dist = "centroids.dist")
tune_feats$choice.keepX

res = tune_feats$choice.keepX # set the optimal values of features to retain
res

keepX = res

final.diablo.model = block.splsda(X = train_data, Y = Y_train, ncomp = 2, 
                                  design = design1, keepX = keepX)


final.diablo.model$weights
final.diablo.model$names

# see selected variables
reacs1 = selectVar(final.diablo.model, block = 'fluxomics', comp = 1)$fluxomics$name
reacs2 = selectVar(final.diablo.model, block = 'fluxomics', comp = 2)$fluxomics$name

mets1 = selectVar(final.diablo.model, block = 'metabolomics', comp = 1)$metabolomics$name
mets2 = selectVar(final.diablo.model, block = 'metabolomics', comp = 2)$metabolomics$name

genes1 = selectVar(final.diablo.model, block = 'RNASeq', comp = 1)$RNASeq$name
genes2 = selectVar(final.diablo.model, block = 'RNASeq', comp = 2)$RNASeq$name

write.csv(reacs2, 
          file="featselect_all_reacs2.csv")

Y = final.diablo.model$Y

# plotting
plotDiablo(final.diablo.model, ncomp = 1)
plotDiablo(final.diablo.model, ncomp = 2)

# PLS plot mas com todos
plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 18), cex = c(2,2, 2), 
        col = c('darkorchid', 'brown1', 'cyan'), cutoff = 0.5)


# correlations between features of the dataset  
corr = circosPlot(final.diablo.model, cutoff = 0.90, line = TRUE,
                  color.blocks= c('cyan3', 'chartreuse1', 'brown1'),
                  color.cor = c("chocolate3","grey20"), size.labels = 1.5, 
                  size.variables = 0.40)

# save corr between features of the dataset
write.csv(corr, 
          file="correlation_all_featselect_NEW.csv")

plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median', ndisplay = 15, size.legend = 1, size.name = 1)
plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'mean', ndisplay = 15, size.legend = 1, size.name = 1)


loads = final.diablo.model$loadings
loads

loads = final.diablo.model$loadings$metabolomics
x = sort(abs(loads[,'comp2']))
names(x)


train_data$metabolomics$shikimic_acid
median(train_data$metabolomics$ferulic_acid)


# model evaluation
perf.diablo = perf(final.diablo.model, validation = 'Mfold', 
                   M = 10, nrepeat = 10, 
                   auc = TRUE) 

# perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate
perf.diablo$error.rate

perf.diablo$AveragedPredict.error.rate
#perf.diablo$WeightedPredict.error.rate

perf.diablo$auc


auc.splsda = auroc(final.diablo.model, roc.block = "RNASeq", 
                   roc.comp = 1, print = FALSE)

auc.splsda = auroc(final.diablo.model, roc.block = "metabolomics", 
                   roc.comp = 1, print = FALSE)

auroc(final.diablo.model, roc.block = "fluxomics", 
      roc.comp = 1, print = FALSE)



# check performance TEST DATA

data_test = list(RNASeq = Xtest_dna,
                 metabolomics = Xtest_met,
                 fluxomics = Xtest_flux)

predict.diablo = predict(final.diablo.model, newdata = data_test)


confusion.mat = get.confusion_matrix(truth = Y_test,
                                     predicted = predict.diablo$WeightedVote$mahalanobis.dist[,1])

get.BER(confusion.mat)

Y_test
