// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include <math.h>
#include <cmath>
#include <stdint.h>

// Based on http://www.csre.iitb.ac.in/~avikb/GNR401/DIP/DIP_401_lecture_7.pdf
// likelihood = (1 / (a * b)) * c
// a = (2pi)^(L/2)
// L: number of classes
// b = det(covarianceMatrix)^(1/2)
// c = exp(d_t * inv(covarianceMatrix) * d)
// d = (X - mu)
// X: feature vector
// mu: mean of feature vector

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List cpp_MLC(const arma::mat X, const arma::vec Y, const CharacterVector levels) {
    const int L = levels.length();
    const double a = std::pow(2.0 * M_PI, L/2.0);
    const double EPSILON = std::pow(DBL_EPSILON, 2/3);
    

    NumericVector k(L);
    auto mu = List::create();
    auto inverseCovarianceMatrices = List::create();

    
    // For each i group
    for (int i = 1; i <= L; i++) {
        // Get Xi submatrix
        arma::uvec ids = arma::find(Y == i);
        arma::mat Xi = X.rows(ids);

        // mu_i
        mu.push_back((arma::mat)arma::mean(Xi));

        // covariance matrix
        arma::mat covarianceMatrix = arma::cov(Xi);

        // Add EPSILON to diagonal of cov to avoid singularity
        covarianceMatrix.diag() += EPSILON; 

        // inverse matrix 
        arma::mat inverseCovarianceMatrix = arma::inv_sympd(covarianceMatrix); 
        inverseCovarianceMatrices.push_back(inverseCovarianceMatrix);

        // determinant of covariance matrix: log_det is more precise than raw determinant
        arma::cx_double logDeterminant = arma::log_det(covarianceMatrix); 
        double det = exp(real(logDeterminant));

        // Constant fraction for each group 1/((2*pi)^(L/2) * determinant(covarianceMatrix_i))
        k[i-1] = 1.0 / (a * std::sqrt(det));
    }    

    auto resultList = List::create(
        _["groups"] = levels,
        _["mu"] = mu,
        _["k"] = k,
        _["inverseCovarianceMatrices"] = inverseCovarianceMatrices
    );
    
    return(resultList);
}


// k * exp(-dist %*% inverseCovarianceMatrix %*% dist_t)
// [[Rcpp::export]]
arma::mat cpp_predict(S4 model, arma::mat X) {
    CharacterVector levels = model.slot("groups");
    const int L = levels.length();    
    const int n = X.n_rows;
    
    List mu = model.slot("mu");
    NumericVector k = model.slot("k");
    List inverseCovarianceMatrices = model.slot("inverseCovarianceMatrices");
    
    arma::mat classLikelihoods(n, L);

    for (int i = 0; i < L; i++) {
        // Retrieve objects from model
        const double k_i = k[i];
        arma::rowvec mu_i = mu[i];
        arma::mat inverseCovarianceMatrix_i = inverseCovarianceMatrices[i];

        // (x - mu) * inverseCovarianceMatrix * transpose(x - mu)
        // transpose work only if classifying feature by feature
        // instead we are using sumByRow(((x - mu) * inverseCovarianceMatrix) % (x-mu))
        // then each row will contain the exponential part of likelihood
        //
        // where % is the element-wise multiplication
        arma::mat distance = X.each_row() - mu_i;
        arma::mat likelihood = arma::sum((distance * inverseCovarianceMatrix_i) % distance, 1);
        likelihood *= -0.5;
        likelihood = arma::exp(likelihood);
        likelihood *= k_i;

        classLikelihoods.col(i) = likelihood;
    }
    
    return(classLikelihoods);
}
