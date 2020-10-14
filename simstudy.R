rm(list=ls())
source("../../mmsca/mmscaScripts/makeData.R")
source("../../mmsca/mmscaScripts/tucker.R")

library(Rcpp)
library(doParallel)

###################################################################################################################
# Simulation study for Tucker and correctly classified

# Note: in this simulation study I forgot to save the raw data, I retrieve this later in a second simulation study
# Later on in this script


# Start simulation study

library(doParallel)
                               

# Explanation: This function performs a "mini" simulation study, given the number of replications (reps)
# And the conditions
simulation_study <- function(reps, n, p, ncomp, sparsity, varianceOfComps, error) {

    datObjectRes <- vector("list", length = reps)
    estimationCCregPCA <- vector("list", length = reps)
    estimationSparsePCA <- vector("list", length = reps)
    tuckerSpca <- vector("list", length = reps)
    tuckerCCregPCA <- vector("list", length = reps)
    corClassSpca <- vector("list", length = reps)
    corClassCCregPCA <- vector("list", length = reps)


    for (i in 1:reps) {
        try({
        # Generate data
        comdis <- matrix(1, p, ncomp)
        comdisspar <- sparsify(comdis, sparsity)
        variances <- makeVariance(varianceOfComps, p = p, error = error)
        datObject <- makeDat(n = n, p = p, ncomp = ncomp, comdisspar, variances = variances)

        resSpca <- spca(datObject$X, K = ncomp, para = apply(comdisspar, 2, sum), type = c("predictor"),
                  sparse = "varnum", use.corr = FALSE, lambda = 1e-6, max.iter = 100000,
                  trace = TRUE, eps.conv = 10^-3)

        resCCRegPCA <- pcaccregCpp(datObject$X, Q = ncomp, apply(comdisspar, 2, function(x){sum(x==0)}), 10000000,
                                   Wstart = matrix(0, ncol(datObject$X), 3), nStarts = 1, tol = 10^-8, printLoss = TRUE)


        # Save results
        estimationSparsePCA[[i]] <- resSpca
        estimationCCregPCA[[i]] <- resCCRegPCA
        tuckerSpca[[i]] <- tuckerCongruence(resSpca$loadings, datObject$P[, 1:ncomp])
        tuckerCCregPCA[[i]] <- tuckerCongruence(resCCRegPCA$W, datObject$P[, 1:ncomp])
        corClassSpca[[i]] <- corClass(datObject$P[, 1:ncomp], resSpca$loadings, comdisspar = comdisspar)
        corClassCCregPCA[[i]] <- corClass(datObject$P[, 1:ncomp], resCCRegPCA$W, comdisspar = comdisspar)
        })
    }

    #NOTE: I forgot to save the simulalted data in datObjectRes, I added that part here later!
    return(list(datObjectRes = datObjectRes, # I added this part after the fact
                estimationCCregPCA = estimationCCregPCA,
                estimationSparsePCA = estimationSparsePCA,
                tuckerSpca = tuckerSpca,
                tuckerCCregPCA = tuckerCCregPCA,
                corClassSpca = corClassSpca,
                corClassCCregPCA = corClassCCregPCA))
}


# Create conditions
nCond <- list(25, 50, 100)
pCond <- list(50)
ncompCond <- list(3)
sparsityCond <- list(c(0.3, 0.3, 0.3), c(0.8,0.8,0.8))
variancesCond  <- list(c(31, 30, 29))
errorCond <- list(0.05, 0.2)

# Create conditions
conditions <- list(nCond, pCond, ncompCond, sparsityCond, variancesCond, errorCond)

conditionCombs <- expand.grid(1:length(nCond), 1:length(pCond), 1:length(ncompCond),
                             1:length(sparsityCond), 1:length(variancesCond), 1:length(errorCond))

nrow(conditionCombs)

conditionCombs


# Create duplicates of the conditions so you can run them multiple times
# So the cores get used more optimally
# NOTE: I do this because: some conditions always take longer than others, so I duplicate them "x" times, 
# and make the number of replicated data sets small, so one core is never overloaded with a task that takes too long

reps <- 10 
x <- 5 

# In total reps*x per conditions
conditionCombsReplicated <- conditionCombs
for (i in 1:(x-1)) {
    conditionCombsReplicated <- rbind(conditionCombsReplicated, conditionCombs)
}

conditionCombsReplicated
(1:nrow(conditionCombsReplicated) %% nrow(conditionCombs)) #these groups belong together

distribute_conditions <- function(i, reps, conditionCombs, conditions) {
               
    set.seed(i) 
    n <- conditions[[1]][[conditionCombs[i, 1]]] 
    p <- conditions[[2]][[conditionCombs[i, 2]]] 
    ncomp <- conditions[[3]][[conditionCombs[i, 3]]] 
    sparsity <- conditions[[4]][[conditionCombs[i, 4]]] 
    varianceOfComps <- conditions[[5]][[conditionCombs[i, 5]]] 
    error <- conditions[[6]][[conditionCombs[i, 6]]] 


    out <- simulation_study(reps = reps, n = n, p = p, ncomp = ncomp, sparsity = sparsity,
                            varianceOfComps = varianceOfComps, error = error)


    saveRDS(out, file=paste("./sim_res/results_", i, sep = "")) 
    return(out) 
}



# Run the simulation study 

########### Start the simulation study
# intialise cores


cores <- detectCores() - 2
cores
cl <- makeCluster(cores, outfile = './sim_res/simstudy_log') # I save output to a log
registerDoParallel(cl)

rm(simulation_study_results)



#Start measuring time
startTime <- Sys.time()

# Run the simulation study

# Note: I use the elasticnet package and my own CCregPCA package, my own package needs to be installed
# Because the code needs to be compiled before destributing over cores
simulation_study_results <- foreach(i = 1:nrow(conditionCombsReplicated), 
                                    .packages = c('elasticnet', 'CCregPCA'), 
                                    .errorhandling = 'pass') %dopar% { 
    source("../../mmsca/mmscaScripts/makeData.R", local = TRUE)
    source("../../mmsca/mmscaScripts/tucker.R", local = TRUE)

    distribute_conditions(i, reps = reps, conditionCombsReplicated, conditions)
}

stopCluster(cl)
endTime <- Sys.time()
timeTaken  <- endTime - startTime
timeTaken


test <- readRDS("./sim_res/results_1") 
test


###################################################################################################################
# Extract only the simulated data

# Note: In this simulation study I extract the data I forgot to save earlier in the previous simulation study
# That is the only purpose of this simulation study

simulation_study2 <- function(reps, n, p, ncomp, sparsity, varianceOfComps, error) {

    datObjectRes <- vector("list", length = reps)
    estimationCCregPCA <- vector("list", length = reps)
    estimationSparsePCA <- vector("list", length = reps)
    tuckerSpca <- vector("list", length = reps)
    tuckerCCregPCA <- vector("list", length = reps)
    corClassSpca <- vector("list", length = reps)
    corClassCCregPCA <- vector("list", length = reps)


    for (i in 1:reps) {
        try({
        comdis <- matrix(1, p, ncomp)
        comdisspar <- sparsify(comdis, sparsity)
        variances <- makeVariance(varianceOfComps, p = p, error = error)
        datObject <- makeDat(n = n, p = p, ncomp = ncomp, comdisspar, variances = variances)

        #resSpca <- spca(datObject$X, K = ncomp, para = apply(comdisspar, 2, sum), type = c("predictor"),
        #          sparse = "varnum", use.corr = FALSE, lambda = 1e-6, max.iter = 100000,
        #          trace = TRUE, eps.conv = 10^-3)

        #resCCRegPCA <- pcaccregCpp(datObject$X, Q = ncomp, apply(comdisspar, 2, function(x){sum(x==0)}), 10000000,
        #                           Wstart = matrix(0, ncol(datObject$X), 3), nStarts = 1, tol = 10^-8, printLoss = TRUE)


        #estimationSparsePCA[[i]] <- resSpca
        #estimationCCregPCA[[i]] <- resCCRegPCA
        #tuckerSpca[[i]] <- tuckerCongruence(resSpca$loadings, datObject$P[, 1:ncomp])
        #tuckerCCregPCA[[i]] <- tuckerCongruence(resCCRegPCA$W, datObject$P[, 1:ncomp])
        #corClassSpca[[i]] <- corClass(datObject$P[, 1:ncomp], resSpca$loadings, comdisspar = comdisspar)
        #corClassCCregPCA[[i]] <- corClass(datObject$P[, 1:ncomp], resCCRegPCA$W, comdisspar = comdisspar)

        datObjectRes[[i]] <- datObject
        })
    }

    return(list(datObjectRes = datObjectRes,
                estimationCCregPCA = estimationCCregPCA,
                estimationSparsePCA = estimationSparsePCA,
                tuckerSpca = tuckerSpca,
                tuckerCCregPCA = tuckerCCregPCA,
                corClassSpca = corClassSpca,
                corClassCCregPCA = corClassCCregPCA))
}


#create conditions
nCond <- list(25, 50, 100)
pCond <- list(50)
ncompCond <- list(3)
sparsityCond <- list(c(0.3, 0.3, 0.3), c(0.8,0.8,0.8))
variancesCond  <- list(c(31, 30, 29))
errorCond <- list(0.05, 0.2)

#create conditions
conditions <- list(nCond, pCond, ncompCond, sparsityCond, variancesCond, errorCond)

conditionCombs <- expand.grid(1:length(nCond), 1:length(pCond), 1:length(ncompCond),
                             1:length(sparsityCond), 1:length(variancesCond), 1:length(errorCond))

nrow(conditionCombs)

conditionCombs


#create duplicates of the conditions so you can run them multiple times
#so the cores get used more optimally
reps <- 10 
x <- 5 

#in total reps*x per conditions
conditionCombsReplicated <- conditionCombs
for (i in 1:(x-1)) {
    conditionCombsReplicated <- rbind(conditionCombsReplicated, conditionCombs)
}

conditionCombsReplicated
(1:nrow(conditionCombsReplicated) %% nrow(conditionCombs)) #these groups belong together

distribute_conditions <- function(i, reps, conditionCombs, conditions) {
               
    set.seed(i) 
    n <- conditions[[1]][[conditionCombs[i, 1]]] 
    p <- conditions[[2]][[conditionCombs[i, 2]]] 
    ncomp <- conditions[[3]][[conditionCombs[i, 3]]] 
    sparsity <- conditions[[4]][[conditionCombs[i, 4]]] 
    varianceOfComps <- conditions[[5]][[conditionCombs[i, 5]]] 
    error <- conditions[[6]][[conditionCombs[i, 6]]] 


    out <- simulation_study2(reps = reps, n = n, p = p, ncomp = ncomp, sparsity = sparsity,
                            varianceOfComps = varianceOfComps, error = error)


    saveRDS(out, file=paste("./sim_res_rawdata/results_", i, sep = "")) 
    return(out) 
}

########### Start the simulation study
# intialise cores


cores <- detectCores() - 2
cores
cl <- makeCluster(cores, outfile = './sim_res/simstudy_log')
registerDoParallel(cl)

rm(simulation_study_results)

simulation_study_results <- foreach(i = 1:nrow(conditionCombsReplicated), 
                                    .packages = c('elasticnet', 'CCregPCA'), 
                                    .errorhandling = 'pass') %dopar% { 
    source("../../mmsca/mmscaScripts/makeData.R", local = TRUE)
    source("../../mmsca/mmscaScripts/tucker.R", local = TRUE)

    distribute_conditions(i, reps = reps, conditionCombsReplicated, conditions)
}

stopCluster(cl)
endTime <- Sys.time()
timeTaken  <- endTime - startTime
timeTaken



####################################################################################################
# Simulation study: BIAS

##### Functions to run the simulation study

simulation_study_bias <- function(reps, n, ncomp, datObject, comdisspar){
    p <- ncol(datObject$Sigma) 

    estimationCCregPCA <- vector("list", length = reps)
    estimationSparsePCA <- vector("list", length = reps)

    for (i in 1:reps) {
        try({

        X <- mvrnorm(n, mu = rep(0, p), datObject$Sigma, empirical=F)

        resSpca <- spca(X, K = ncomp, para = apply(comdisspar, 2, sum), type = c("predictor"),
                  sparse = "varnum", use.corr = FALSE, lambda = 1e-6, max.iter = 100000,
                  trace = TRUE, eps.conv = 10^-3)

        resCCRegPCA <- pcaccregCpp(X, Q = ncomp, apply(comdisspar, 2, function(x){sum(x==0)}), 10000000,
                                   Wstart = matrix(0, ncol(datObject$X), 3), nStarts = 1, tol = 10^-8, printLoss = TRUE)


        estimationSparsePCA[[i]] <- resSpca
        estimationCCregPCA[[i]] <- resCCRegPCA

        })
    }

    return(list(datObjectRes = datObject,
                estimationCCregPCA = estimationCCregPCA,
                estimationSparsePCA = estimationSparsePCA))
}


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
        

    set.seed(i+4) # Always the same data objects
    comdis <- matrix(1, p, ncomp)
    comdisspar <- sparsify(comdis, sparsity)
    variances <- makeVariance(varianceOfComps, p = p, error = error)
    datObject <- makeDat(n = n, p = p, ncomp = ncomp, comdisspar, variances = variances)
    dataObjectList[[i]] <- datObject
    comdissparList[[i]] <- comdisspar
}


# Remove not used conditions from the conditions list
conditions[[6]] <- NULL 
conditions[[5]] <- NULL
conditions[[4]] <- NULL
conditions[[2]] <- NULL

# Attach dataObjectList and comdissparList to the conditions
conditions[[3]] <- dataObjectList
conditions[[4]] <- comdissparList

# Remove all unused conditions from conditionCombs
conditionCombs <- conditionCombs[, -1*c(6,5,4,2)]
# Add the condition numbers for the data object and sparse structure
conditionCombs <- cbind(conditionCombs, 1:nrow(conditionCombs), 1:nrow(conditionCombs))

conditionCombs


distribute_conditions_bias <- function(i, reps, conditionCombs, conditions) {
               
    set.seed(i) 
    n <- conditions[[1]][[conditionCombs[i, 1]]] 
    ncomp <- conditions[[2]][[conditionCombs[i, 2]]] 
    datObject <- conditions[[3]][[conditionCombs[i, 3]]] 
    comdisspar <- conditions[[4]][[conditionCombs[i, 4]]] 

    out <- simulation_study_bias(reps = reps, n = n, ncomp = ncomp, datObject = datObject,
                                 comdisspar = comdisspar)

    saveRDS(out, file=paste("./sim_res/results_bias_", i, sep = "")) 
    return(out) 
}

# Create duplicates of the conditions so you can run them multiple times
# so the cores get used more optimally
reps <- 50 
x <- 100 

# In total reps*x per conditions
conditionCombsReplicated <- conditionCombs
for (i in 1:(x-1)) {
    conditionCombsReplicated <- rbind(conditionCombsReplicated, conditionCombs)
}

conditionCombsReplicated
(1:nrow(conditionCombsReplicated) %% nrow(conditionCombs)) #these groups belong together



########################################################################################

########### Start the simulation study
# Intialise cores

cores <- detectCores() - 1
cores
cl <- makeCluster(cores, outfile = './sim_res/simstudy_log')
registerDoParallel(cl)


# Start measuring time
startTime <- Sys.time()

# Run the simulation study
simulation_study_results <- foreach(i = 1:nrow(conditionCombsReplicated), 
                                    .packages = c('elasticnet', 'CCregPCA'), 
                                    .errorhandling = 'pass') %dopar% { 

    source("../../mmsca/mmscaScripts/makeData.R", local = TRUE)
    source("../../mmsca/mmscaScripts/tucker.R", local = TRUE)

    distribute_conditions_bias(i, reps = reps, conditionCombsReplicated, conditions)
}


stopCluster(cl)
endTime <- Sys.time()
timeTaken  <- endTime - startTime
timeTaken

check <- readRDS("./sim_res/results_bias_700")
names(check)
length(check$estimationCCregPCA)
check

