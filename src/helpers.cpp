// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include <math.h>
#include <cmath>
#include <stdint.h>

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
List cpp_MLC(const arma::mat X, const arma::vec Y, const CharacterVector levels) {
    const int L = levels.length();
    const double a = std::pow(2.0 * M_PI, L/2.0);
    const double EPSILON = std::pow(DBL_EPSILON, 2/3);
    

    NumericVector k(L);
    auto mu = List::create();
    auto invCovs = List::create();

    
    for (int i = 1; i <= L; i++) {
        arma::uvec ids = arma::find(Y == i); // Find indices
        arma::mat subX = X.rows(ids);
        mu.push_back((arma::mat)arma::mean(subX));


        arma::mat covMat = arma::cov(subX);
        covMat.diag() += EPSILON;
        arma::mat invCov = arma::inv_sympd(covMat);
        invCovs.push_back(invCov);

        // Determinants
        arma::cx_double logDet = arma::log_det(covMat);
        double det = exp(real(logDet));
        double k_i = 1.0 / (a * std::sqrt(det));
        k[i-1] = k_i;
        

    }    

    auto resultList = List::create(
        _["groups"] = levels,
        _["mu"] = mu,
        _["k"] = k,
        _["invCovs"] = invCovs
    );
    
    return(resultList);
}


// k * exp(-dist %*% invCov %*% dist_t)
// [[Rcpp::export]]
arma::mat cpp_predict(S4 model, arma::mat X) {
    CharacterVector levels = model.slot("groups");
    const int L = levels.length();    
    const int n = X.n_rows;
    
    List mu = model.slot("mu");
    NumericVector k = model.slot("k");
    List invCovs = model.slot("invCovs");
    
    arma::mat classLikelihoods(n, L);

    for (int i = 0; i < L; i++) {
        const double k_i = k[i];
        arma::rowvec mu_i = mu[i];
        arma::mat invCov_i = invCovs[i];

        arma::mat distance = X.each_row() - mu_i;
        arma::mat likelihood = arma::sum((distance * invCov_i) % distance,1);
        likelihood *= -0.5;
        likelihood = arma::exp(likelihood);
        likelihood.transform( [k_i](double val) { return (val * k_i); } );

        classLikelihoods.col(i) = likelihood;
    }
    
    return(classLikelihoods);
}