##############################################################
# Create the bias variance MSE table from the paper
##############################################################
rm(list=ls())

##############################################################
# Load the data

# Script to make the bias tables
test <- readRDS("./sim_res/results_bias_1")
names(test)
test$estimationCCregPCA
test$tuckerSpca
test$corClassCCregPCA
test$corClassSpca

# Create conditions for the bias sim study
nCond <- list(25, 100)
pCond <- list(50)
ncompCond <- list(3)
sparsityCond <- list(c(0.3, 0.3, 0.3), c(0.8,0.8,0.8))
variancesCond  <- list(c(30, 30, 30))
errorCond <- list(0.05, 0.2)

conditions <- list(nCond, pCond, ncompCond, sparsityCond, variancesCond, errorCond)

conditionCombs <- expand.grid(1:length(nCond), 1:length(pCond), 1:length(ncompCond),
                             1:length(sparsityCond), 1:length(variancesCond), 1:length(errorCond))

dataObjectList <- vector("list", length = nrow(conditionCombs))
comdissparList <- vector("list", length = nrow(conditionCombs))
source("../../mmsca/mmscaScripts/makeData.R")
for (i in 1:nrow(conditionCombs)) {
    
    n <- conditions[[1]][[conditionCombs[i, 1]]] 
    p <- conditions[[2]][[conditionCombs[i, 2]]] 
    ncomp <- conditions[[3]][[conditionCombs[i, 3]]] 
    sparsity <- conditions[[4]][[conditionCombs[i, 4]]] 
    varianceOfComps <- conditions[[5]][[conditionCombs[i, 5]]] 
    error <- conditions[[6]][[conditionCombs[i, 6]]] 
        

    set.seed(i+4) #always the same data objects
    comdis <- matrix(1, p, ncomp)
    comdisspar <- sparsify(comdis, sparsity)
    variances <- makeVariance(varianceOfComps, p = p, error = error)
    datObject <- makeDat(n = n, p = p, ncomp = ncomp, comdisspar, variances = variances)
    dataObjectList[[i]] <- datObject
    comdissparList[[i]] <- comdisspar
}


# Create duplicates of the conditions so you can run them multiple times
# So the cores get used more optimally
reps <- 50 
x <- 100 

# In total reps*x per conditions
conditionCombsReplicated <- conditionCombs
for (i in 1:(x-1)) {
    conditionCombsReplicated <- rbind(conditionCombsReplicated, conditionCombs)
}


conditionCombsReplicated
grouping <- (0:(nrow(conditionCombsReplicated)-1) %% nrow(conditionCombs)) + 1 #these groups belong together
grouping


# Combine the estimates from all the files
resListCCregPCA <- vector("list", length = nrow(conditionCombs))
resListSparsePCA <- vector("list", length = nrow(conditionCombs))

for (a in 1:nrow(conditionCombs)) {
    res <- matrix(NA, reps, 12)
    files <- 1:nrow(conditionCombsReplicated)  
    files <- files[grouping == a]
    print(files)

    for (i in 1:x) {
        fileName <- paste("./sim_res/results_bias_", files[i], sep = "")
        loadedFile <- readRDS(fileName)

        resListCCregPCA[[a]]  <- c(resListCCregPCA[[a]], loadedFile$estimationCCregPCA)
        resListSparsePCA[[a]] <- c(resListSparsePCA[[a]], loadedFile$estimationSparsePCA)
    }
}

#################################################################################################################
# Calculate the bias


# Function to determine the order of What in relation to W
# use the tucker congruence to mix and match

# Important note: this is the function I used to rearrange What given W
# I took into account the sign in an absolute wrong way
# I created rearrangeWhat2 function that should fix this issue, I compared results
# where I used rearrangeWhat and rearrangeWhat2, and there are NO difference in results! 
# In comming simulation studies rearrangeWhat2 should be used (the next user please verify this function)

rearrangeWhat <- function(What, W) {

    rearrangedWhat <- matrix(NA, nrow(What), ncol(What))
    combinationList <- combinat::permn(1:(ncol(What)))
    candidates <- vector("list", length = length(combinationList))
    fixWhat <- What

    for (i in 1:length(combinationList)) {
        What <- fixWhat[, combinationList[[i]]]
        tucCon <- rep(NA, ncol(What))
        for (j in 1:ncol(What)) {
            tucCon[j] <- What[, j] %*% W[, j]  / sqrt(sum(What[, j]^2) * sum(W[, j]^2))
        }
        candidates[[i]] <- tucCon 
    }
    best <- lapply(candidates, function(x){mean(abs(x))})
    best <- which.max(best)
    What <- fixWhat[, combinationList[[best]] ] %*% diag(sign(candidates[[best]]))
    return(What)
}

source("./tucker.R")
rearrangeWhat2 <- function(What, W) {

    rearrangedWhat <- matrix(NA, nrow(What), ncol(What))
    combinationList <- combinat::permn(1:(ncol(What)))
    candidates <- vector("list", length = length(combinationList))
    candidatesSign <- vector("list", length = length(combinationList))
    changeInSign <- c(1, -1)
    placeHolder <- rep(NA, 2)
    fixWhat <- What

    for (i in 1:length(combinationList)) {
        What <- fixWhat[, combinationList[[i]]]
        tucCon <- rep(NA, ncol(What))
        signage <- rep(NA, ncol(What))
        for (j in 1:ncol(What)) {
            for (a in c(1, 2)) {
                placeHolder[a] <- changeInSign[a]*What[, j] %*% W[, j]  / sqrt(sum(What[, j]^2) * sum(W[, j]^2))
            }
            tucCon[j] <- placeHolder[which.max(placeHolder)]
            signage[j] <- changeInSign[which.max(placeHolder)]
        }
        candidates[[i]] <- tucCon 
        candidatesSign[[i]] <- signage 
    }
    best <- lapply(candidates, function(x){mean(abs(x))})
    best <- which.max(best)
    What <- fixWhat[, combinationList[[best]] ] %*% diag(candidatesSign[[best]])
    return(What)
}


# Function to estimate the MAB ccreg
estimate_bias_ccreg <- function(results, W) {

    reps <- length(results)
    Q <- ncol(results[[1]]$W)
    W <- W[, 1:Q]
    sumWhatCCregPCA <- matrix(0, nrow(W), ncol(W))

    for (i in 1:reps) {
        W_CCregPCA <- rearrangeWhat(results[[i]]$W, W)
        sumWhatCCregPCA  <- sumWhatCCregPCA + W_CCregPCA 
    }
    meanWhatCCregPCA <- sumWhatCCregPCA / reps
    sumWhatCCregPCA <- matrix(0, nrow(W), ncol(W))
    for (i in 1:reps) {
        W_CCregPCA <- rearrangeWhat(results[[i]]$W, W)
        sumWhatCCregPCA  <- sumWhatCCregPCA + (W_CCregPCA - meanWhatCCregPCA)^2
    }
    varmat <- sumWhatCCregPCA / (reps - 1)
    MAB_CCregPCA <- sum(abs(meanWhatCCregPCA - W)) / (ncol(W) * nrow(W))

    return(list(MAB = MAB_CCregPCA, varmat = varmat))

}

# Function to estimate the MAB sparsepca
estimate_bias_spca <- function(results, W) {

    reps <- length(results)
    Q <- ncol(results[[1]]$loadings)
    W <- W[, 1:Q]
    sumWhatSpca <- matrix(0, nrow(W), ncol(W))

    for (i in 1:reps) {
        W_Spca <- rearrangeWhat(results[[i]]$loadings, W)
        sumWhatSpca <- sumWhatSpca + W_Spca
    }
    meanWhatSpca <-  sumWhatSpca / reps
    sumWhatSpca <- matrix(0, nrow(W), ncol(W))
    for (i in 1:reps) {
        W_Spca <- rearrangeWhat(results[[i]]$loadings, W)
        sumWhatSpca  <- sumWhatSpca + (W_Spca - meanWhatSpca)^2
    }
    varmat <- sumWhatSpca / (reps - 1)
    MAB_Spca <- sum(abs(meanWhatSpca - W)) / (ncol(W) * nrow(W))
    return(list(MAB = MAB_Spca, varmat = varmat))
}

##################################

# Function to estimate the MSE ccreg
estimate_mse_ccreg <- function(results, W) {

    reps <- length(results)
    Q <- ncol(results[[1]]$W)
    W <- W[, 1:Q]

    sumWhatCCregPCA <- matrix(0, nrow(W), ncol(W))
    for (i in 1:reps) {
        W_CCregPCA <- rearrangeWhat(results[[i]]$W, W)
        sumWhatCCregPCA  <- sumWhatCCregPCA + (W_CCregPCA - W)^2
    }
    varmat <- sumWhatCCregPCA / (reps - 1)
    MSE <- sum(varmat) / (ncol(W) * nrow(W))

    return(list(MSE = MSE))
}

# Function to estimate the MSE spca
estimate_mse_spca <- function(results, W) {

    reps <- length(results)
    Q <- ncol(results[[1]]$loadings)
    W <- W[, 1:Q]

    sumWhatSpca <- matrix(0, nrow(W), ncol(W))
    for (i in 1:reps) {
        W_Spca <- rearrangeWhat(results[[i]]$loadings, W)
        sumWhatSpca  <- sumWhatSpca + (W_Spca - W)^2
    }
    varmat <- sumWhatSpca / (reps - 1)
    MSE <- sum(varmat) / (ncol(W) * nrow(W))

    return(list(MSE = MSE))
}


# Estimate bias ccreg
MAB_CCregPCA <- c(
estimate_bias_ccreg(resListCCregPCA[[1]], dataObjectList[[1]]$P)$MAB,
estimate_bias_ccreg(resListCCregPCA[[2]], dataObjectList[[2]]$P)$MAB,
estimate_bias_ccreg(resListCCregPCA[[3]], dataObjectList[[3]]$P)$MAB,
estimate_bias_ccreg(resListCCregPCA[[4]], dataObjectList[[4]]$P)$MAB,
estimate_bias_ccreg(resListCCregPCA[[5]], dataObjectList[[5]]$P)$MAB,
estimate_bias_ccreg(resListCCregPCA[[6]], dataObjectList[[6]]$P)$MAB,
estimate_bias_ccreg(resListCCregPCA[[7]], dataObjectList[[7]]$P)$MAB,
estimate_bias_ccreg(resListCCregPCA[[8]], dataObjectList[[8]]$P)$MAB
)

# Estimate variance ccreg
mean_var_CCregPCA <- c(
mean(estimate_bias_ccreg(resListCCregPCA[[1]], dataObjectList[[1]]$P)$varmat),
mean(estimate_bias_ccreg(resListCCregPCA[[2]], dataObjectList[[2]]$P)$varmat),
mean(estimate_bias_ccreg(resListCCregPCA[[3]], dataObjectList[[3]]$P)$varmat),
mean(estimate_bias_ccreg(resListCCregPCA[[4]], dataObjectList[[4]]$P)$varmat),
mean(estimate_bias_ccreg(resListCCregPCA[[5]], dataObjectList[[5]]$P)$varmat),
mean(estimate_bias_ccreg(resListCCregPCA[[6]], dataObjectList[[6]]$P)$varmat),
mean(estimate_bias_ccreg(resListCCregPCA[[7]], dataObjectList[[7]]$P)$varmat),
mean(estimate_bias_ccreg(resListCCregPCA[[8]], dataObjectList[[8]]$P)$varmat)
)

# Estimate mse ccreg
mean_mse_CCregPCA <- c(
estimate_mse_ccreg(resListCCregPCA[[1]], dataObjectList[[1]]$P)$MSE,
estimate_mse_ccreg(resListCCregPCA[[2]], dataObjectList[[2]]$P)$MSE,
estimate_mse_ccreg(resListCCregPCA[[3]], dataObjectList[[3]]$P)$MSE,
estimate_mse_ccreg(resListCCregPCA[[4]], dataObjectList[[4]]$P)$MSE,
estimate_mse_ccreg(resListCCregPCA[[5]], dataObjectList[[5]]$P)$MSE,
estimate_mse_ccreg(resListCCregPCA[[6]], dataObjectList[[6]]$P)$MSE,
estimate_mse_ccreg(resListCCregPCA[[7]], dataObjectList[[7]]$P)$MSE,
estimate_mse_ccreg(resListCCregPCA[[8]], dataObjectList[[8]]$P)$MSE
)

# Estimate bias spca
MAB_Spca <- c(
estimate_bias_spca(resListSparsePCA[[1]], dataObjectList[[1]]$P)$MAB,
estimate_bias_spca(resListSparsePCA[[2]], dataObjectList[[2]]$P)$MAB,
estimate_bias_spca(resListSparsePCA[[3]], dataObjectList[[3]]$P)$MAB,
estimate_bias_spca(resListSparsePCA[[4]], dataObjectList[[4]]$P)$MAB,
estimate_bias_spca(resListSparsePCA[[5]], dataObjectList[[5]]$P)$MAB,
estimate_bias_spca(resListSparsePCA[[6]], dataObjectList[[6]]$P)$MAB,
estimate_bias_spca(resListSparsePCA[[7]], dataObjectList[[7]]$P)$MAB,
estimate_bias_spca(resListSparsePCA[[8]], dataObjectList[[8]]$P)$MAB
)

# Estimate variance spca
mean_var_Spca <- c(
mean(estimate_bias_spca(resListSparsePCA[[1]], dataObjectList[[1]]$P)$varmat),
mean(estimate_bias_spca(resListSparsePCA[[2]], dataObjectList[[2]]$P)$varmat),
mean(estimate_bias_spca(resListSparsePCA[[3]], dataObjectList[[3]]$P)$varmat),
mean(estimate_bias_spca(resListSparsePCA[[4]], dataObjectList[[4]]$P)$varmat),
mean(estimate_bias_spca(resListSparsePCA[[5]], dataObjectList[[5]]$P)$varmat),
mean(estimate_bias_spca(resListSparsePCA[[6]], dataObjectList[[6]]$P)$varmat),
mean(estimate_bias_spca(resListSparsePCA[[7]], dataObjectList[[7]]$P)$varmat),
mean(estimate_bias_spca(resListSparsePCA[[8]], dataObjectList[[8]]$P)$varmat)
)

# Estimate mse spca
mean_mse_Spca <- c(
estimate_mse_spca(resListSparsePCA[[1]], dataObjectList[[1]]$P)$MSE,
estimate_mse_spca(resListSparsePCA[[2]], dataObjectList[[2]]$P)$MSE,
estimate_mse_spca(resListSparsePCA[[3]], dataObjectList[[3]]$P)$MSE,
estimate_mse_spca(resListSparsePCA[[4]], dataObjectList[[4]]$P)$MSE,
estimate_mse_spca(resListSparsePCA[[5]], dataObjectList[[5]]$P)$MSE,
estimate_mse_spca(resListSparsePCA[[6]], dataObjectList[[6]]$P)$MSE,
estimate_mse_spca(resListSparsePCA[[7]], dataObjectList[[7]]$P)$MSE,
estimate_mse_spca(resListSparsePCA[[8]], dataObjectList[[8]]$P)$MSE
)


nCond <- list(25, 100)
pCond <- list(50)
ncompCond <- list(3)
sparsityCond <- list(c(0.3, 0.3, 0.3), c(0.8,0.8,0.8))
variancesCond  <- list(c(30, 30, 30))
errorCond <- list(0.05, 0.2)

conditions <- list(nCond, pCond, ncompCond, sparsityCond, variancesCond, errorCond)

conditionCombs <- expand.grid(1:length(nCond), 1:length(pCond), 1:length(ncompCond),
                             1:length(sparsityCond), 1:length(variancesCond), 1:length(errorCond))

bias_var_mse <- as.data.frame(cbind(c(MAB_CCregPCA,
                    MAB_Spca),
                    c(mean_var_CCregPCA,
                    mean_var_Spca),
                    c(mean_mse_CCregPCA,
                    mean_mse_Spca)))

bias_var_mse <- cbind(bias_var_mse, c(rep("CCregPCA", nrow(conditionCombs)), rep("Spca", nrow(conditionCombs))))

bias_var_mse <- cbind(bias_var_mse, rbind(conditionCombs, conditionCombs))
colnames(bias_var_mse) <- c("MAB", "mVar", "mMSE", "Method", "I", "J", "Q", "Sparsity", "Variance", "Error")

bias_var_mse$Method <- as.factor(bias_var_mse$Method)
bias_var_mse$Sparsity <- as.factor(bias_var_mse$Sparsity)
bias_var_mse$Variance <- as.factor(bias_var_mse$Variance)
bias_var_mse$Error <- as.factor(bias_var_mse$Error)
bias_var_mse$I <- as.factor(bias_var_mse$I) 

levels(bias_var_mse$Sparsity) <- c("Sparsity 30%", "Sparsity 80%")
levels(bias_var_mse$Variance) <- c("Variances 30 30 30")
levels(bias_var_mse$Error) <- c("Error 5%", "Error 20%")
levels(bias_var_mse$I) <- c("25", "100")

library(xtable)

# This produces the latex skeleton table!
# Values need to be copy and pasted in!!
xtable(bias_var_mse, digits = 5)
tab <- ftable(bias_var_mse$Method, bias_var_mse$Error, bias_var_mse$I, bias_var_mse$Sparsity, row.vars = 1)
xtableFtable(tab, align = rep("c", 9))
xtableFtable(tab)

# Produce values to be copy and pasted into the table.
bias_var_mse <-  bias_var_mse[order(bias_var_mse$Sparsity), ]
bias_var_mse <-  bias_var_mse[order(bias_var_mse$I), ]
bias_var_mse <-  bias_var_mse[order(bias_var_mse$Error), ]

noquote(paste(" & ", sprintf("%.4f", round(bias_var_mse$MAB[1:10 %% 2 != 0], 4)), sep = ""))
noquote(paste(" & ", sprintf("%.4f", round(bias_var_mse$MAB[1:10 %% 2 == 0], 4)), sep = ""))

noquote(paste(" & ", sprintf("%.4f", round(bias_var_mse$mVar[1:10 %% 2 != 0], 4)), sep = ""))
noquote(paste(" & ", sprintf("%.4f", round(bias_var_mse$mVar[1:10 %% 2 == 0], 4)), sep = ""))

noquote(paste(" & ", sprintf("%.4f", round(bias_var_mse$mMSE[1:10 %% 2 != 0], 4)), sep = ""))
noquote(paste(" & ", sprintf("%.4f", round(bias_var_mse$mMSE[1:10 %% 2 == 0], 4)), sep = ""))


