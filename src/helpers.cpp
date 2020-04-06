// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include <math.h>
#include <cmath>

// Based on http://www.csre.iitb.ac.in/~avikb/GNR401/DIP/DIP_401_lecture_7.pdf
// likelihood = (1 / (a * b)) * c
// a = (2pi)^(L/2)
// L: number of classes
// b = det(covMatrix)^(1/2)
// c = exp(d_t * inv(covMatrix) * d)
// d = (X - mu)
// X: feature vector
// mu: mean of feature vector

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function

using namespace Rcpp;

// available from R
//
// [[Rcpp::export]]
List cpp_MLC(arma::mat X, arma::rowvec Y, CharacterVector levels) {
    int L = levels.length();

    NumericVector k(L);
    auto mu = List::create();
    auto invCovs = List::create();

    double a = std::pow(2.0 * M_PI, L/2.0);
    Rprintf("a = %lf\n", a);
    
    for (int i = 1; i <= L; i++) {
        arma::uvec ids = find(Y == i); // Find indices
        arma::mat subX = X.rows(ids);
        mu.push_back((arma::mat)arma::mean(subX));

        arma::mat covMat = arma::cov(subX);
        double b = std::sqrt(arma::det(covMat));
        Rprintf("i:%d  b:%lf\n", i, b);
        k[i-1] = 1.0 / (a*b);
        invCovs.push_back(arma::inv_sympd(covMat));
    }    

    auto resultList = List::create(
        _["groups"] = levels,
        _["mu"] = mu,
        _["k"] = k,
        _["invCovs"] = invCovs
    );
    
    return(resultList);
}

