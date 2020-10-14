##################################################################
# Function for the Tucker congruence coefficient
##################################################################

#check the tucker congruence between vector a and all columns of B
#returns the highest tucker congruence between a and the column of be
#also returns the location of the column in B
tuckerCongruenceOnevsMany <- function(a, B){
    res <- -Inf
    changeInSign <- c(1, -1)
    placeHolder <- rep(NA, 2)
    tucCon <- rep(NA, ncol(B))
    for( i in 1:ncol(B) ){
        for(j in c(1, 2)){
            placeHolder[j] <- changeInSign[j]*a %*% B[, i]  / sqrt(sum(a^2) * sum(B[, i]^2))
            }
        tucCon[i] <- max(placeHolder)
    }
    return(list(tucker=max(tucCon), where=which.max(tucCon)))
}

A <- matrix(1:16, 4, 4)
B <- matrix(rnorm(16), 4, 4)

combinationList <- combinat::permn(1:(ncol(A)))
res <- rep(NA, length(combinationList))

for (i in combinationList) {
    print(i)
}



#Average tucker congruence coefficient
#of all columns of A and B
tuckerCongruence <- function(A, B){                                                                 
    combinationList <- combinat::permn(1:(ncol(A)))
    res <- -Inf
    fixA <- A
    changeInSign <- c(1, -1)
    placeHolder <- rep(NA, 2)

        for(combination in combinationList){
            A <- fixA[, combination]
            tucCon <- rep(NA, ncol(A))
               for( i in 1:ncol(A) ){
                    for(j in c(1, 2)){
                    placeHolder[j] <- changeInSign[j]*A[, i] %*% B[, i]  / sqrt(sum(A[, i]^2) * sum(B[, i]^2))
                }
                tucCon[i] <- max(placeHolder)
            }
            candidate <- mean(tucCon)
            if(!is.nan(candidate) && candidate > res){
                res <- candidate
            }
        }
    return(res)
}


countZero <- function(A){
  return(sum( A == 0 ) / length(A))
}

correctlyClassified <- function(A, B){
  counter <- 0
  for( i in 1:nrow(A) ){
    for( j in 1:ncol(A) ){
      if(A[i, j] == B[i, j]){
        counter <- counter + 1
      }
      if(A[i, j] != 0 && B[i, j] != 0){
        counter <- counter + 1
      }
    }
  }
  return(counter / (nrow(A)*ncol(A)))
}


#function that checks where the distinctive components are
#then assesses whether all coefficients from the zero blocks,
#are actually zero
#the common components are indicated with NA
disIdentified <- function(W, What, comdis){
    ncomp <- ncol(What)
    distinctiveIdentified <- rep(NA, ncomp)
    for(i in 1:ncomp){
        where <- tuckerCongruenceOnevsMany(W[, i], What)$where
        coefsThatShouldBeZero <- which(comdis[, i] == 0)
        if(length(coefsThatShouldBeZero) != 0) {
            distinctiveIdentified[i] <- sum(What[coefsThatShouldBeZero, where] == 0) / length(coefsThatShouldBeZero)
        } else {
            distinctiveIdentified[i] <- NA 
        }
    }
    return(distinctiveIdentified)
}


#function that finds whether zero's vs non-zero's are correctly classified
corClass <- function(W, What, comdisspar){
    ncomp <- ncol(What)

    zeroCoefFound <- rep(NA, ncomp)
    nonZeroCoefFound <- rep(NA, ncomp)
    for(i in 1:ncomp){
        where <- tuckerCongruenceOnevsMany(W[, i], What)$where
        coefsThatShouldBeZero <- which(comdisspar[, i] == 0)
        coefsThatShouldBeNonZero <- which(comdisspar[, i] == 1)

        zeroCoefFound[i] <- sum(What[coefsThatShouldBeZero, where] == 0) / length(coefsThatShouldBeZero)
        nonZeroCoefFound[i] <- sum(What[coefsThatShouldBeNonZero, where] != 0) / length(coefsThatShouldBeNonZero)

    }
    return(list(zeroCoefFound = zeroCoefFound, nonZeroCoefFound = nonZeroCoefFound))
}

#Function that checks where the common components are
#then assesses whether at least one variable in each block is selected
#the distinctive components are indicated with NA
comIdentified <- function(W, What, comdis, groups){
    ncomp <- ncol(What)
    ngroups <- length(groups)
    cumsumgroups <- c(0, cumsum(groups))
    commonIdentified <- rep(NA, ncomp)
    for (i in 1:ncomp) {
        where <- tuckerCongruenceOnevsMany(W[, i], What)$where
        coefsThatShouldBeNonZero <- which(comdis[, i] == 1)
        if (length(coefsThatShouldBeNonZero) == nrow(comdis)) {
            atLeastOneNonZeroCoef <- 0
            for (j in 1:ngroups) {
                indices <- (cumsumgroups[j] + 1):cumsumgroups[j+1]
                atLeastOneNonZeroCoef <- atLeastOneNonZeroCoef + ifelse(sum(What[indices, where]) != 0, 1, 0)
            }
            commonIdentified[i] <- atLeastOneNonZeroCoef / ngroups
        } else {
            commonIdentified[i] <- NA 
        }
    }
    return(commonIdentified)
}

