#include <RcppArmadillo.h>
#include <cmath>


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void updatewCpp(arma::mat& W, const arma::mat& P, const int& Q, const arma::mat& XtX, 
    const arma::vec& nzeroes, const double& a) {

    arma::mat Sxy = XtX * P;

    for (int i = 0; i < Q; i++) {
        W.col(i) = W.col(i) - (1/a) * (XtX * W.col(i) - Sxy.col(i)); 
        //add cardinality constraints
        arma::uvec indices = arma::sort_index(arma::abs(W.col(i)));
        for (int j = 0; j < nzeroes[i]; j++) {
            W(indices[j], i) = 0;
        }
    }

}

// [[Rcpp::export]]
Rcpp::List pcaccregCpp(const arma::mat& X, const int& Q, const arma::vec& nzeroes, 
        const int &itr, arma::mat Wstart, int nStarts, double tol, bool printLoss) {


    const arma::mat XtX = X.t() * X;
    arma::mat U, V, U2, V2, P, W;
    arma::vec D, D2;
    Rcpp::List ret;

    double prevLoss, curLoss, minLossGlobal = arma::datum::inf;
    bool converged;

    arma::svd(U, D, V, X);
    double a = std::pow(arma::max(D), 2);

    for (int k = 0; k < nStarts; k++) {
        prevLoss = arma::datum::inf;
        curLoss = std::pow(10, 100);
        converged = false;

        /* First run will be warm, or with the custom starting values */
        if(arma::accu(Wstart) == 0){
            arma::svd(U, D, V, X);
            W = V.cols(0, Q - 1);
            W += arma::randu<arma::mat>(size(W)) * k;
        } else {
            /* If the starting value matrix does not sum to zero,
             * meaning the user wants custom starting values, do:
             */
            W = Wstart;
            W += arma::randu<arma::mat>(size(W)) * k;
        }

        for (int i = 0; i < itr; i++) {
            Rcpp::checkUserInterrupt();
            //check condition if true break else continue
            if (prevLoss - curLoss < tol) {
                Rcpp::Rcout << "converged" << "\n";
                converged = true;
                break; 
            }
            //print loss 
            if (printLoss && i % 1000 == 0) {
                Rcpp::Rcout << "Start: " << k+1 << " at itr: " << i << " at loss val: " << curLoss << "\n";
            }

            //procruste rotation least squares P given W
            arma::svd(U2, D2, V2, XtX * W);
            //arma::svds(U2, D2, V2, arma::sp_mat(XtX * W), Q);
            //P = U2 * V2.t(); 
            P = U2.cols(0, Q-1) * V2.t(); 
            
            //updateW
            updatewCpp(W, P, Q, XtX, nzeroes, a);
            
            prevLoss = curLoss;
            curLoss =  arma::accu(arma::pow((X - X * W * P.t()), 2));
        }

        if(curLoss < minLossGlobal){
            minLossGlobal = curLoss;

            //set small elements to zero
            ret["W"] = W;
            ret["P"] = P; 
            ret["loss"] = minLossGlobal;
            ret["converged"] = converged;
        } 
    }

    return ret;
}


