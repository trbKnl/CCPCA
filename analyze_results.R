#########################################################################
# Script to make the plots in the paper

rm(list=ls())

#########################################################################
# Read in the data

#some testing 
test <- readRDS("./sim_res/results_1")

test$tuckerCCregPCA
test$tuckerSpca
test$corClassCCregPCA
test$corClassSpca


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
grouping <- (0:(nrow(conditionCombsReplicated)-1) %% nrow(conditionCombs)) + 1 #these groups belong together

results_to_analyse <- vector("list", length = nrow(conditionCombs))

for (a in 1:nrow(conditionCombs)) {
    res <- matrix(NA, reps, 12)
    files <- 1:nrow(conditionCombsReplicated)  
    files <- files[grouping == a]
    print(files)
    resList <- vector("list", length = x)

    for (i in 1:x) {
        #fileName <- paste("./simulation_study_modelselection_results/model_selection_result_", files[i], sep = "")
        fileName <- paste("./sim_res/results_", files[i], sep = "")
        loadedFile <- readRDS(fileName)

        for (j in 1:reps) {

            res[j, 1] <- loadedFile$tuckerCCregPCA[[j]]
            res[j, 2] <- loadedFile$tuckerSpca[[j]]

            res[j, 3] <- mean(loadedFile$corClassCCregPCA[[j]]$zeroCoefFound)
            res[j, 4] <- mean(loadedFile$corClassSpca[[j]]$zeroCoefFound)

            res[j, 5] <- mean(loadedFile$corClassCCregPCA[[j]]$nonZeroCoefFound)
            res[j, 6] <- mean(loadedFile$corClassSpca[[j]]$nonZeroCoefFound)

            res[j, 7] <- conditionCombs[a, 1]
            res[j, 8] <- conditionCombs[a, 2]
            res[j, 9] <- conditionCombs[a, 3]
            res[j, 10] <- conditionCombs[a, 4]
            res[j, 11] <- conditionCombs[a, 5]
            res[j, 12] <- conditionCombs[a, 6]

        }
        resList[[i]] <- res
    }
    resultsForCondition <-  do.call(rbind, resList)
    results_to_analyse[[a]] <- cbind(resultsForCondition, a)
}

results_to_analyse


results_to_analyse <- do.call(rbind, results_to_analyse)
results_to_analyse <- as.data.frame(results_to_analyse)
results_to_analyse 



colnames(results_to_analyse)  <- c("tucker_CCregPCA", "tucker_Spca", "zeroCoefFound_CCregPCA", "zeroCoefFound_Spca", 
                                   "nonZeroCoefFound_CCregPCA", "nonZeroCoefFound_Spca", "n", "p", "ncomp", "sparsity", 
                                   "variances", "error", "conditionNumber")

results_to_analyse <- reshape(results_to_analyse, varying = 1:6, sep = "_", timevar = "method", direction = "long")

results_to_analyse$error <- as.factor(results_to_analyse$error)
results_to_analyse$n <- as.factor(results_to_analyse$n)
results_to_analyse$sparsity <- as.factor(results_to_analyse$sparsity)

levels(results_to_analyse$sparsity) <- c("Sparsity 30%", "Sparsity 80%")
levels(results_to_analyse$error) <- c("Error 5%", "Error 20%")
levels(results_to_analyse$n) <- c("25", "50", "100")

results_to_analyse$id <- NULL
results_to_analyse$totalCoefFound <- NA

TF <- results_to_analyse$sparsity == "Sparsity 30%"
results_to_analyse$totalCoefFound[TF] <-  (results_to_analyse$zeroCoefFound[TF] * (150*0.3) + (150*0.7) * results_to_analyse$nonZeroCoefFound[TF]) / 150

TF <- results_to_analyse$sparsity == "Sparsity 80%"
results_to_analyse$totalCoefFound[TF] <-  (results_to_analyse$zeroCoefFound[TF] * (150*0.8) + (150*0.2) * results_to_analyse$nonZeroCoefFound[TF]) / 150

library(ggplot2)

# scripts to make the boxplots
names(results_to_analyse)

# rename CCregPCA to CCPCA 
results_to_analyse$method[results_to_analyse$method == "CCregPCA"]  <- "CCPCA"
results_to_analyse$method[results_to_analyse$method == "Spca"]  <- "sparse PCA"

############################### tucker plot
plot <- ggplot(results_to_analyse, aes(y = tucker, x = method))
plot + geom_boxplot(aes(fill = n), outlier.size=0.3, lwd=0.2) +
    facet_grid(rows = vars(error), cols = vars(sparsity)) +
    geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
    guides(fill=guide_legend(title="I")) +
    xlab("Method") + 
    scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
    ylab("Tucker congruence") +
    coord_fixed(ratio = 26/5) +
    theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=10),
                axis.text.y = element_text(size=10),
                strip.text = element_text(size=10),
                legend.text = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(0.5, "cm"))

#ggsave("~/work/pca_ccreg/main/plots/Tucker.pdf", height = 10, width = 10, unit='in', dpi = 300)

#crop whitespace with pdfcrop
#system("pdfcrop --margin 5  ~/work/pca_ccreg/main/plots/Tucker.pdf ~/work/pca_ccreg/main/plots/Tucker.pdf")

##############################################################################
##############################################################################
##############################################################################
##############################################################################
# weights found plot
results_to_analyse2 <- results_to_analyse 
results_to_analyse2$id  <- NULL


rownames(results_to_analyse2) <- NULL
colnames(results_to_analyse2)[10] <- "percentage_zeroCoefFound" 
colnames(results_to_analyse2)[11] <- "percentage_nonZeroCoefFound"
colnames(results_to_analyse2)[12] <- "percentage_totalCoefFound"
head(results_to_analyse2)

results_to_analyse2$id <- 1:1200
head(results_to_analyse2)

results_to_analyse2 <- reshape(results_to_analyse2, sep ="_", direction = "long", varying = 10:12, timevar="found") 
results_to_analyse2
results_to_analyse2$found

results_to_analyse2$found <- factor(results_to_analyse2$found, levels = c("nonZeroCoefFound",  "zeroCoefFound", "totalCoefFound"), ordered = TRUE)
levels(results_to_analyse2$found)
levels(results_to_analyse2$found) <- c("Non-zero weights" , "Zero weights", "All weights")
rownames(results_to_analyse2) <-  NULL

plot <- ggplot(results_to_analyse2, aes(y = percentage, x = method ))
plot + geom_boxplot(aes(fill = found), outlier.size =0.3, lwd=0.2,
                   position = "dodge2") +
    facet_grid(error + sparsity ~ n,
              scales="free_x",
              space="free_x",
              drop=TRUE) +
    guides(fill=guide_legend(title="Weights")) +
    xlab("Method") + 
    ylab("Proportion of correctly identified weights") +
    scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3", "#D3D3D3")) +
    theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=10),
                axis.text.y = element_text(size=10),
                strip.text = element_text(size=10),
                legend.text = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(0.5, "cm"))



#ggsave("~/work/pca_ccreg/main/plots/WeightsFound.pdf", height = 10, width = 10, unit='in', dpi = 300)
#
##crop whitespace with pdfcrop
#system("pdfcrop --margin 5  ~/work/pca_ccreg/main/plots/WeightsFound.pdf ~/work/pca_ccreg/main/plots/WeightsFound.pdf")
#

summary(results_to_analyse2)
head(results_to_analyse2)


#########################################################################
# Add correctly estimating weight based on chance to the boxplots

# 15 is 30% of 50
n <- 50 # Total number of weights
k <- 15 # Total number of non-zero  weights
result2 <- matrix(NA, 2, k+1)


for (i in 0:k) {
    # Percentage total correct
    result2[1, i+1] <- (n - (k-i + k-i)) / n 

    # Probablity of that occuring
    result2[2, i+1] <-  (choose(n-k, k-i) * choose(k, i)) /choose(n, k) 
}

result2[2, ] <- round(result2[2, ] * 10000)

out <- c()
for (i in 1:ncol(result2)) {
    out <- c(out, rep(result2[1, i], result2[2, i])    )
}

len <- length(out)

n <- rep("Chance", len)
p <- rep(1, len) 
ncomp <- rep(1, len) 
sparsity <- rep("Sparsity 30%", len)
variances <- rep(1, len)
error5 <- rep("Error 5%", len)
error20 <- rep("Error 20%", len)
conditionNumber <- rep(1, len)
method <- rep("Chance", len)
tucker <- rep(1, len)
id <- rep(1, len)
found <- rep("All weights by chance", len)

add1 <- data.frame(n, p, ncomp, sparsity, variances, error = error5, conditionNumber, method, tucker, id, found, percentage = out)
add2 <- data.frame(n, p, ncomp, sparsity, variances, error = error20, conditionNumber, method, tucker, id, found, percentage = out)


results_to_analyse2 <- rbind(results_to_analyse2, add1, add2)

# 40 is 80% of 50
n <- 50
k <- 40 
result2 <- matrix(NA, 2, k+1)

# NOTE: this contains values that cannot occur, for example only 1 correct CANNOT happen
# but can safely be ignored because they have probablity zero of occuring
for (i in 0:k) {
    result2[1, i+1] <- (n - (k-i + k-i)) / n
    result2[2, i+1] <-  (choose(n-k, k-i) * choose(k,i)) /choose(n, k)
}


result2[2, ] <- round(result2[2, ] * 10000)
result2
out <- c()
for (i in 1:ncol(result2)) {
    out <- c(out, rep(result2[1, i], result2[2, i])    )
}
len <- length(out)

n <- rep("Chance", len)
p <- rep(1, len) 
ncomp <- rep(1, len) 
sparsity <- rep("Sparsity 80%", len)
variances <- rep(1, len)
error5 <- rep("Error 5%", len)
error20 <- rep("Error 20%", len)
conditionNumber <- rep(1, len)
method <- rep("Chance", len)
tucker <- rep(1, len)
id <- rep(1, len)
found <- rep("All weights by chance", len)

add1 <- data.frame(n, p, ncomp, sparsity, variances, error = error5, conditionNumber, method, tucker, id, found, percentage = out)
add2 <- data.frame(n, p, ncomp, sparsity, variances, error = error20, conditionNumber, method, tucker, id, found, percentage = out)

results_to_analyse2 <- rbind(results_to_analyse2, add1, add2)



plot <- ggplot(results_to_analyse2, aes(y = percentage, x = method ))
p <- plot + geom_boxplot(aes(fill = found), outlier.size =0.3, lwd=0.2,
                position = position_dodge2(preserve = "single")) +
    facet_grid(error + sparsity ~ n,
              scales="free_x",
              space="free_x",
              drop=TRUE) +
    guides(fill=guide_legend(title="Weights")) +
    xlab("Method") + 
    ylab("Proportion of correctly identified weights") +
    scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3", "#ffffff")) +
    theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=10),
                axis.text.y = element_text(size=10),
                strip.text = element_text(size=10),
                legend.text = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(0.5, "cm"))

p

#ggsave("~/work/pca_ccreg/main/plots/WeightsFound.pdf", plot = p, height = 8, width = 8, unit='in', dpi = 300)

#crop whitespace with pdfcrop
#system("pdfcrop --margin 5  ~/work/pca_ccreg/main/plots/WeightsFound.pdf ~/work/pca_ccreg/main/plots/WeightsFound.pdf")


#####################################################################
# calculate average absolute value of simulate component weight
# Note: this number is used in the paper

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
grouping <- (0:(nrow(conditionCombsReplicated)-1) %% nrow(conditionCombs)) + 1 #these groups belong together

results_to_analyse <- vector("list", length = nrow(conditionCombs))

test <- readRDS("./sim_res_rawdata/results_1")

for (a in 1:nrow(conditionCombs)) {
    res <- matrix(NA, reps, 12)
    files <- 1:nrow(conditionCombsReplicated)  
    files <- files[grouping == a]
    print(files)
    resList <- vector("list", length = x)

    for (i in 1:x) {
        #fileName <- paste("./simulation_study_modelselection_results/model_selection_result_", files[i], sep = "")
        fileName <- paste("./sim_res_rawdata/results_", files[i], sep = "")
        loadedFile <- readRDS(fileName)

        for (j in 1:reps) {

            res[j, 1] <- mean(abs(test$datObjectRes[[j]]$P[1:3, ]))

        }
        resList[[i]] <- res
    }
    resultsForCondition <-  do.call(rbind, resList)
    results_to_analyse[[a]] <- cbind(resultsForCondition, a)
}

results_to_analyse <- do.call(rbind, results_to_analyse)
mean(results_to_analyse[, 1])
