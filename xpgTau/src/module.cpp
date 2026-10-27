// #include <iostream>
// #include <cmath>
// #include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ==========================
// Autoreg.Rates.vec (C++)
// ==========================
// [[Rcpp::export]]
arma::vec Autoreg_Rates_vec_cpp(int r, const arma::mat& xs, const arma::vec& thetas, double G) {
  int n = xs.n_rows; // the length of rat is the row number of xs
  arma::vec rat(n, fill::zeros);

  switch (r) {
  case 1: rat = thetas(0) * (G - xs.col(3)) % xs.col(2); break;       // DNA+P2->DNA.P2; % means element-wise multiplication
  case 2: rat = thetas(1) * xs.col(3); break;                          // DNA.P2->DNA+P2
  case 3: rat = thetas(2) * (G - xs.col(3)); break;                    // DNA->DNA+RNA
  case 4: rat = thetas(3) * xs.col(0); break;                          // RNA->RNA+P
  case 5: rat = thetas(4) * xs.col(1) % (xs.col(1) - 1.0) / 2.0; break;// 2P->P2
  case 6: rat = thetas(5) * xs.col(2); break;                          // P2->2P
  case 7: rat = thetas(6) * xs.col(0); break;                          // RNA->0
  case 8: rat = thetas(7) * xs.col(1); break;                          // P->0
  default: break;
  }
  return rat;
}

// [[Rcpp::export]]
arma::vec SEIR_Rates_coupled_vec(
    int r,
    const arma::mat& xs,
    const arma::vec& thetas,
    const arma::mat& coupled_A,
    const arma::vec& Ns
) {
  int D = xs.n_cols / 4;
  int M1 = xs.n_rows;
  
  arma::uvec index_S = regspace<uvec>(0, D - 1);
  arma::uvec index_E = regspace<uvec>(D, 2 * D - 1);
  arma::uvec index_I = regspace<uvec>(2 * D, 3 * D - 1);
  
  // Compute rates
  arma::mat R1 = thetas(0) * xs.cols(index_S) % (xs.cols(index_I) * coupled_A.t()); // (M+1) x D
  arma::mat R2 = thetas(1) * xs.cols(index_E);
  arma::mat R3 = thetas(2) * xs.cols(index_I);
  
  arma::vec rat(M1, arma::fill::zeros);
  
  for (int d = 0; d < D; ++d) {
    if (r == d + 1) {                        // infection event
      rat = R1.col(d) / Ns(d);
      break;
    } else if (r == D + d + 1) {             // exposed → infectious
      rat = R2.col(d);
      break;
    } else if (r == 2 * D + d + 1) {         // infectious → removed
      rat = R3.col(d);
      break;
    }
  }
  
  return rat;
}


// [[Rcpp::export]]
arma::vec SIR_Rates_coupled_vec(
    int r,
    const arma::mat& xs,
    const arma::vec& thetas,
    const arma::mat& coupled_A,
    const arma::vec& Ns
) {
  int D = xs.n_cols / 3;
  int M1 = xs.n_rows;

  arma::uvec index_S = regspace<uvec>(0, D - 1);
  arma::uvec index_I = regspace<uvec>(D, 2 * D - 1);

  // thetas(0) = beta
  // thetas(1) = gamma

  // Infection rates:
  // beta * S_d * sum_j coupled_A(d,j) * I_j / N_d
  arma::mat R1 = thetas(0) * xs.cols(index_S) % (xs.cols(index_I) * coupled_A.t());

  // Removal/recovery rates:
  // gamma * I_d
  arma::mat R2 = thetas(1) * xs.cols(index_I);

  arma::vec rat(M1, arma::fill::zeros);

  for (int d = 0; d < D; ++d) {
    if (r == d + 1) {
      // S_d -> I_d
      rat = R1.col(d) / Ns(d);
      break;
    } else if (r == D + d + 1) {
      // I_d -> R_d
      rat = R2.col(d);
      break;
    }
  }
  
  return rat;
}


// [[Rcpp::export]]
arma::vec SIS_Rates_independent_vec(
    int r,
    const arma::mat& xs,
    const arma::vec& thetas,
    const arma::vec& Ns
) {
  // xs: (M+1) x D
  // columns are I_1, ..., I_D
  //
  // thetas(0) = beta
  // thetas(1) = gamma

  int D = xs.n_cols;
  int M1 = xs.n_rows;

  arma::vec rat(M1, arma::fill::zeros);

  for (int d = 0; d < D; ++d) {

    if (r == d + 1) {

      // Infection:
      // I_d -> I_d + 1
      //
      // rate = beta * (N_d - I_d) * I_d / N_d

      arma::vec I_d = xs.col(d);
      arma::vec S_d = Ns(d) - I_d;

      rat = thetas(0) * (S_d % I_d) / Ns(d);

      break;

    } else if (r == D + d + 1) {

      // Recovery:
      // I_d -> I_d - 1
      //
      // rate = gamma * I_d

      rat = thetas(1) * xs.col(d);

      break;
    }
  }

  return rat;
}


// [[Rcpp::export]]
IntegerVector make_sequence(int start, int end) {
  if (start <= end) return Rcpp::seq(start, end);
  else return Rcpp::rev(Rcpp::seq(end, start));
}


// [[Rcpp::export]]
arma::vec update_proposals_cpp(const IntegerVector& indices, const std::string& direction,
                                const IntegerVector& sorted_indices,
                                const arma::vec& Mu,
                                arma::vec tilde_N_proposals,
                                double delta) {

  int len = indices.size();
  if (direction == "forward") {
    for (int j = 0; j < len; ++j) {
      int current = sorted_indices[indices[j] - 1]; // 1-based to 0-based
      int next_i = sorted_indices[indices[j] + 1 -1]; // 1-based to 0-based
      double lambda = (Mu[next_i] - Mu[current]) * (1 - delta);
      tilde_N_proposals[next_i] = tilde_N_proposals[current] + R::rpois(lambda);
    }
  } else { // backward
    for (int j = 0; j < len; ++j) {
      int current = sorted_indices[indices[j] - 1];
      int prev_i = sorted_indices[indices[j] - 1 - 1];
      if (Mu[current] > 0) {
        tilde_N_proposals[prev_i] = R::rbinom(tilde_N_proposals[current], Mu[prev_i] / Mu[current]);
      } else {
        tilde_N_proposals[prev_i] = 0;
      }
    }
  }
  return tilde_N_proposals;
}

// ==========================
// propagate (C++)
// ==========================
// [[Rcpp::export]]
List propagate_rcpp(int l,
                   arma::mat X_ancestors,
                   const arma::mat& X_ref_local, // dim(X_ref_local) = (u-1)*(time_steps[l+1]-time_steps[l])
                   const arma::mat& N_ref_local, // dim(N_ref_local) = v*(time_steps[l+1]-time_steps[l])
                   double delta,
                   int M,
                   const arma::vec& time_steps,
                   const arma::mat& A,
                   double G,
                   double tau,
                   const arma::vec& thetas,
                   Rcpp::Function rate_func) {

  int u = A.n_rows;
  int v = A.n_cols;
  int dt = time_steps(l) - time_steps(l - 1);

  // Allocate 3D cubes
  arma::cube N_proposals_local(M + 1, v, dt, arma::fill::value(arma::datum::nan));
  arma::cube X_proposals_local(M + 1, u - 1, dt, arma::fill::value(arma::datum::nan)); 

  // Initialize first row with references
    // ---- Equivalent to N.proposals_local[1, , ] <- N.ref_local ----
  for (int k = 0; k < dt; ++k) {
    // N_ref_local: v x dt  -> N_proposals_local.slice(k).row(0) should match length v
    N_proposals_local.slice(k).row(0) = trans(N_ref_local.col(k));
    // X_ref_local: (u-1) x dt
    X_proposals_local.slice(k).row(0) = trans(X_ref_local.col(k));
  }
  
  // Main loop
  for (int k = 0; k < dt; ++k) {
    // tilde.N.ref
    arma::vec tilde_N_ref(v);
    for (int j = 0; j < v; ++j){
      tilde_N_ref(j) = R::rbinom(N_ref_local(j, k), 1 - delta);
      } // generate tilde_N for each reaction for k

    // tilde_N_proposals
    arma::mat tilde_N_proposals(M+1, v, fill::value(datum::nan));
    tilde_N_proposals.row(0) = tilde_N_ref.t(); // fill in the first row of tilde_N_proposals

    
    // update_proposals per column
    for (int i = 0; i < v; ++i) {
      // arma::vec Mu_col = Autoreg_Rates_vec_cpp(i + 1, X_ancestors, thetas, G) * tau; // compute Mu_col for the i-th column
      arma::vec Mu_col = Rcpp::as<arma::vec> (rate_func(i + 1, X_ancestors, thetas, G)) * tau;
      
      Mu_col.replace(datum::nan, 0.0); // replace NA with 0
      arma::vec tilde_N_proposals_col = tilde_N_proposals.col(i); // the i-th col of tilde_N_proposals
      
      arma::uvec sorted_col = arma::sort_index(Mu_col); 
      IntegerVector sorted_indices = seq_len(M + 1);
      sorted_indices = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(Rcpp::IntegerVector(sorted_col.begin(), sorted_col.end()))); // 0-based indices
      
      int c = std::find(sorted_indices.begin(), sorted_indices.end(), 0) - sorted_indices.begin() + 1; // find c for the i-th column; 1-based
    
      if (c == 1) {
        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, M), "forward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
      } else if (c == (M + 1)) {
        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, 2), "backward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
      } else {
        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, 2), "backward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, M), "forward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
      }
      tilde_N_proposals.col(i) = tilde_N_proposals_col;
      for (int row = 1; row < M+1; ++row){
        N_proposals_local(row, i, k) = tilde_N_proposals(row, i) + R::rpois(delta * Mu_col[row]); // fill proposals with Poisson updates
      }
    }

    // X update
    arma::mat N_trunc = N_proposals_local.slice(k);
    arma::mat X_potential(M+1, u, fill::value(datum::nan));

    X_potential.submat(0, 0, 0, u - 2) = X_ref_local.col(k).t();
    arma::mat A_sub = A.rows(0, u - 2);

    for (int row = 1; row < M+1; ++row) {
      arma::rowvec sumv = X_ancestors.row(row) + N_trunc.row(row) * A_sub.t();
      X_potential(row, span(0, u - 2)) = round(sumv);
    }
    X_potential.col(u - 1) = G - X_potential.col(u - 2);

    // Handle negative entries
    for (int row = 1; row < M+1; ++row) {
      if (any(X_potential.row(row) < 0)) {
        X_proposals_local.slice(k).row(row).cols(0, u-2) = X_potential.row(row).cols(0, u-2);
        X_proposals_local.tube(row, 0, row, u - 2).fill(datum::nan);
        N_proposals_local.tube(row, 0, row, v - 1).fill(datum::nan);
      } else {
        X_proposals_local.slice(k).row(row).cols(0, u-2) = X_potential.row(row).cols(0, u-2);
      }
    }

    // update X.ancestors
    X_ancestors = X_proposals_local.slice(k);
  }
  return List::create(
    Named("X.proposals_local") = X_proposals_local,
    Named("N.proposals_local") = N_proposals_local
  );
}


// [[Rcpp::export]]
List propagate_rcpp_SEIR(int l,
                   arma::mat X_ancestors,
                   const arma::mat& X_ref_local, // dim(X_ref_local) = u*(time_steps[l+1]-time_steps[l])
                   const arma::mat& N_ref_local, // dim(N_ref_local) = v*(time_steps[l+1]-time_steps[l])
                   double delta,
                   int M,
                   const arma::vec& time_steps,
                   const arma::mat& A,
                   const arma::vec& Ns,
                   double tau,
                   const arma::vec& thetas,
                   const arma::mat& coupled_A) {

  int u = A.n_rows;
  int v = A.n_cols;
  int dt = time_steps(l) - time_steps(l - 1);

  // Allocate 3D cubes
  arma::cube N_proposals_local(M + 1, v, dt, arma::fill::value(arma::datum::nan));
  arma::cube X_proposals_local(M + 1, u, dt, arma::fill::value(arma::datum::nan)); 

  // Initialize first row with references
    // ---- Equivalent to N.proposals_local[1, , ] <- N.ref_local ----
  for (int k = 0; k < dt; ++k) {
    // N_ref_local: v x dt  -> N_proposals_local.slice(k).row(0) should match length v
    N_proposals_local.slice(k).row(0) = trans(N_ref_local.col(k));
    // X_ref_local: (u-1) x dt
    X_proposals_local.slice(k).row(0) = trans(X_ref_local.col(k));
  }
  
  // Main loop
  for (int k = 0; k < dt; ++k) {
    // tilde.N.ref
    arma::vec tilde_N_ref(v);
    for (int j = 0; j < v; ++j){
      tilde_N_ref(j) = R::rbinom(N_ref_local(j, k), 1 - delta);
      } // generate tilde_N for each reaction for k

    // tilde_N_proposals
    arma::mat tilde_N_proposals(M+1, v, fill::value(datum::nan));
    tilde_N_proposals.row(0) = tilde_N_ref.t(); // fill in the first row of tilde_N_proposals

    
    // update_proposals per column
    for (int i = 0; i < v; ++i) {
      // arma::vec Mu_col = Autoreg_Rates_vec_cpp(i + 1, X_ancestors, thetas, G) * tau; // compute Mu_col for the i-th column
      arma::vec Mu_col = SEIR_Rates_coupled_vec(i + 1, X_ancestors, thetas, coupled_A, Ns) * tau;
      
      
      Mu_col.replace(datum::nan, 0.0); // replace NA with 0
      arma::vec tilde_N_proposals_col = tilde_N_proposals.col(i); // the i-th col of tilde_N_proposals
      
      arma::uvec sorted_col = arma::sort_index(Mu_col); 
      IntegerVector sorted_indices = seq_len(M + 1);
      sorted_indices = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(Rcpp::IntegerVector(sorted_col.begin(), sorted_col.end()))); // 0-based indices
      
      int c = std::find(sorted_indices.begin(), sorted_indices.end(), 0) - sorted_indices.begin() + 1; // find c for the i-th column; 1-based
    
      if (c == 1) {
        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, M), "forward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
      } else if (c == (M + 1)) {
        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, 2), "backward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
      } else {
        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, 2), "backward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, M), "forward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
      }
      tilde_N_proposals.col(i) = tilde_N_proposals_col;
      for (int row = 1; row < M+1; ++row){
        N_proposals_local(row, i, k) = tilde_N_proposals(row, i) + R::rpois(delta * Mu_col[row]); // fill proposals with Poisson updates
      }
    }

    // X update
    arma::mat N_trunc = N_proposals_local.slice(k);
    arma::mat X_potential(M+1, u, fill::value(datum::nan));

    X_potential.row(0) = X_ref_local.col(k).t();
    // arma::mat A_sub = A.rows(0, u - 3);

    for (int row = 1; row < M+1; ++row) {
      arma::rowvec sumv = X_ancestors.row(row) + N_trunc.row(row) * A.t();
      X_potential.row(row) = round(sumv);
    }

    // Handle negative entries
    for (int row = 1; row < M+1; ++row) {
      if (any(X_potential.row(row) < 0)) {
        X_proposals_local.slice(k).row(row).cols(0, u-1) = X_potential.row(row).cols(0, u-1);
        X_proposals_local.tube(row, 0, row, u - 1).fill(datum::nan);
        N_proposals_local.tube(row, 0, row, v - 1).fill(datum::nan);
      } else {
        X_proposals_local.slice(k).row(row).cols(0, u-1) = X_potential.row(row).cols(0, u-1);
      }
    }

    // update X.ancestors
    X_ancestors = X_proposals_local.slice(k);
  }
  return List::create(
    Named("X.proposals_local") = X_proposals_local,
    Named("N.proposals_local") = N_proposals_local
  );
}


// [[Rcpp::export]]
List propagate_rcpp_SIR(
    int l,
    arma::mat X_ancestors,
    const arma::mat& X_ref_local,
    const arma::mat& N_ref_local,
    double delta,
    int M,
    const arma::vec& time_steps,
    const arma::mat& A,
    const arma::vec& Ns,
    double tau,
    const arma::vec& thetas,
    const arma::mat& coupled_A) {
    
    int u = A.n_rows;
    int v = A.n_cols;
    int dt = time_steps(l) - time_steps(l - 1);

    // Allocate 3D cubes
  arma::cube N_proposals_local( M + 1, v, dt, arma::fill::value(arma::datum::nan));

  arma::cube X_proposals_local(M + 1, u, dt, arma::fill::value(arma::datum::nan));


  // ----------------------------------------------------------
  // Initialise reference particle
  // ----------------------------------------------------------
  for (int k = 0; k < dt; ++k) {
  
    // Reference reaction counts
    N_proposals_local.slice(k).row(0) = trans(N_ref_local.col(k));

    // Reference states
    X_proposals_local.slice(k).row(0) = trans(X_ref_local.col(k));
  }

  // ----------------------------------------------------------
  // Main tau-leap loop
  // ----------------------------------------------------------
  for (int k = 0; k < dt; ++k) {

    // --------------------------------------------------------
    // Generate tilde N for reference particle
    // --------------------------------------------------------
    arma::vec tilde_N_ref(v);

    for (int j = 0; j < v; ++j) {

      tilde_N_ref(j) = R::rbinom( N_ref_local(j, k), 1.0 - delta);
    }

    // --------------------------------------------------------
    // Shared/coupled reaction counts
    // --------------------------------------------------------
    arma::mat tilde_N_proposals(M + 1, v, arma::fill::value(arma::datum::nan));

    tilde_N_proposals.row(0) = tilde_N_ref.t();


    // --------------------------------------------------------
    // Update proposals reaction by reaction
    // --------------------------------------------------------
    for (int i = 0; i < v; ++i) {

      // Main difference from SEIR:
      // use the coupled SIR rate function
      arma::vec Mu_col = SIR_Rates_coupled_vec(i + 1, X_ancestors, thetas, coupled_A, Ns) * tau;


      // Replace NA rates by zero
      Mu_col.replace(arma::datum::nan, 0.0);

      arma::vec tilde_N_proposals_col = tilde_N_proposals.col(i);


      // Sort particles according to their Poisson means
      arma::uvec sorted_col = arma::sort_index(Mu_col);

      IntegerVector sorted_indices = seq_len(M + 1);

      sorted_indices = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(Rcpp::IntegerVector(sorted_col.begin(), sorted_col.end())));


      // Position of reference particle in sorted ordering
      // c is 1-based here
      int c = std::find(sorted_indices.begin(), sorted_indices.end(), 0) - sorted_indices.begin() + 1;


      // ------------------------------------------------------
      // Generate shared component
      // ------------------------------------------------------
      if (c == 1) {

        tilde_N_proposals_col = update_proposals_cpp( make_sequence(c, M), "forward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);

      } else if (c == (M + 1)) {

        tilde_N_proposals_col = update_proposals_cpp( make_sequence(c, 2), "backward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);

      } else {

        tilde_N_proposals_col = update_proposals_cpp( make_sequence(c, 2), "backward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);

        tilde_N_proposals_col = update_proposals_cpp( make_sequence(c, M), "forward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
      }

      tilde_N_proposals.col(i) = tilde_N_proposals_col;

      // ------------------------------------------------------
      // Add private Poisson component
      // ------------------------------------------------------
      for (int row = 1; row < M + 1; ++row) {
        N_proposals_local(row, i, k) = tilde_N_proposals(row, i) + R::rpois(delta * Mu_col(row));
      }
    }

    // --------------------------------------------------------
    // Update states
    // --------------------------------------------------------
    arma::mat N_trunc = N_proposals_local.slice(k);

    arma::mat X_potential(M + 1, u, arma::fill::value(arma::datum::nan));


    // Reference particle follows reference path
    X_potential.row(0) = X_ref_local.col(k).t();


    // Other particles
    for (int row = 1; row < M + 1; ++row) {
      arma::rowvec sumv = X_ancestors.row(row) + N_trunc.row(row) * A.t();

      X_potential.row(row) = round(sumv);
    }

    // --------------------------------------------------------
    // Remove invalid particles with negative states
    // --------------------------------------------------------
    for (int row = 1; row < M + 1; ++row) {

      if (any(X_potential.row(row) < 0)) {

        X_proposals_local.slice(k).row(row).cols(0, u - 1) = X_potential.row(row).cols(0, u - 1);

        X_proposals_local.tube(row, 0, row, u - 1).fill(arma::datum::nan);

        N_proposals_local.tube(row, 0, row, v - 1).fill(arma::datum::nan);

      } else {

        X_proposals_local.slice(k).row(row).cols(0, u - 1) = X_potential.row(row).cols(0, u - 1);
      }
    }


    // --------------------------------------------------------
    // Move to next tau step
    // --------------------------------------------------------
    X_ancestors = X_proposals_local.slice(k);
  }

  return List::create(
    Named("X.proposals_local") = X_proposals_local,

    Named("N.proposals_local") = N_proposals_local
  );
}


// [[Rcpp::export]]
List propagate_rcpp_SIS(
    int l,
    arma::mat X_ancestors,
    const arma::mat& X_ref_local,
    const arma::mat& N_ref_local,
    double delta,
    int M,
    const arma::vec& time_steps,
    const arma::mat& A,
    const arma::vec& Ns,
    double tau,
    const arma::vec& thetas
) {

  int u = A.n_rows;   // u = D
  int v = A.n_cols;   // v = 2D
  int dt = time_steps(l) - time_steps(l - 1);


  // ----------------------------------------------------------
  // Allocate proposal arrays
  // ----------------------------------------------------------

  arma::cube N_proposals_local(M + 1, v, dt, arma::fill::value(arma::datum::nan));

  arma::cube X_proposals_local(M + 1, u, dt, arma::fill::value(arma::datum::nan));


  // ----------------------------------------------------------
  // Initialise reference particle
  // ----------------------------------------------------------

  for (int k = 0; k < dt; ++k) {

    // N_ref_local: v x dt
    N_proposals_local.slice(k).row(0) = trans(N_ref_local.col(k));

    // X_ref_local: u x dt
    X_proposals_local.slice(k).row(0) = trans(X_ref_local.col(k));
  }

  // ----------------------------------------------------------
  // Main tau-leap loop
  // ----------------------------------------------------------

  for (int k = 0; k < dt; ++k) {
    // --------------------------------------------------------
    // Shared part of reference reaction counts
    // --------------------------------------------------------

    arma::vec tilde_N_ref(v);

    for (int j = 0; j < v; ++j) {
    tilde_N_ref(j) = R::rbinom( N_ref_local(j, k), 1 - delta);
    }
    
    // --------------------------------------------------------
    // Coupled proposal reaction counts
    // --------------------------------------------------------

    arma::mat tilde_N_proposals(M + 1, v, fill::value(datum::nan));

    tilde_N_proposals.row(0) = tilde_N_ref.t();

    // --------------------------------------------------------
    // Work through reaction channels
    // --------------------------------------------------------

    for (int i = 0; i < v; ++i) {
      // ======================================================
      // SIS rates
      // ======================================================
      arma::vec Mu_col = SIS_Rates_independent_vec(i + 1, X_ancestors, thetas, Ns) * tau;
      
      // Dead/invalid particle
      Mu_col.replace(datum::nan, 0.0);


      arma::vec tilde_N_proposals_col = tilde_N_proposals.col(i);

      // ------------------------------------------------------
      // Sort particles according to mu
      // ------------------------------------------------------

      arma::uvec sorted_col = arma::sort_index(Mu_col);

      IntegerVector sorted_indices = seq_len(M + 1);

      sorted_indices = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(Rcpp::IntegerVector(sorted_col.begin(), sorted_col.end())));


      // Position of reference particle in sorted ordering
      // c is 1-based
      int c = std::find(sorted_indices.begin(), sorted_indices.end(), 0) - sorted_indices.begin() + 1;

      // ------------------------------------------------------
      // Generate shared component
      // ------------------------------------------------------

      if (c == 1) {

        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, M), "forward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);

      } else if (c == (M + 1)) {

        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, 2),"backward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);

      } else {

        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, 2), "backward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);

        tilde_N_proposals_col = update_proposals_cpp(make_sequence(c, M), "forward", sorted_indices, Mu_col, tilde_N_proposals_col, delta);
      }

      tilde_N_proposals.col(i) = tilde_N_proposals_col;

      // ------------------------------------------------------
      // Add private Poisson component
      // ------------------------------------------------------

      for (int row = 1; row < M + 1; ++row) {
        N_proposals_local(row, i, k) = tilde_N_proposals(row, i) + R::rpois(delta * Mu_col(row));
      }
    }
    // --------------------------------------------------------
    // Update SIS states
    // --------------------------------------------------------

    arma::mat N_trunc = N_proposals_local.slice(k);

    arma::mat X_potential(M + 1, u, fill::value(datum::nan));

    // Reference state
    X_potential.row(0) = X_ref_local.col(k).t();

    // Non-reference particles
    for (int row = 1; row < M + 1; ++row) {
    
    arma::rowvec sumv = X_ancestors.row(row) + N_trunc.row(row) * A.t();
    X_potential.row(row) = round(sumv);
    
    }
    // --------------------------------------------------------
    // Check validity
    //
    // For SIS we need
    //
    //        0 <= I_d <= N_d
    //
    // not merely I_d >= 0.
    // --------------------------------------------------------

    for (int row = 1; row < M + 1; ++row) {
      bool invalid = false;
      for (int d = 0; d < u; ++d) {
        if (X_potential(row, d) < 0 || X_potential(row, d) > Ns(d)) {
          invalid = true;
          break;
        }
      }


      if (invalid) {
        X_proposals_local.slice(k).row(row).cols(0, u - 1) = X_potential.row(row).cols(0, u - 1);
        X_proposals_local.tube(row, 0, row, u - 1).fill(datum::nan);
        N_proposals_local.tube(row, 0, row, v - 1).fill(datum::nan);

      } else {
        X_proposals_local.slice(k).row(row).cols(0, u - 1) = X_potential.row(row).cols(0, u - 1);
      }
    }
    // --------------------------------------------------------
    // Move to next tau step
    // --------------------------------------------------------

    X_ancestors = X_proposals_local.slice(k);
  }
  
  return List::create(
    Named("X.proposals_local") = X_proposals_local,
    Named("N.proposals_local") = N_proposals_local
  );
}


// [[Rcpp::export]]
List propagate_rcpp_SEIR_single(
    int l,
    arma::mat X_ancestors,
    const arma::mat& X_ref_local,
    const arma::mat& N_ref_local,
    double delta,
    int M,
    const arma::vec& time_steps,
    const arma::mat& A,
    double N,
    double tau,
    const arma::vec& thetas
) {
  arma::mat coupled_A(1, 1, arma::fill::ones);
  arma::vec Ns(1);
  Ns(0) = N;

  return propagate_rcpp_SEIR(
    l,
    X_ancestors,
    X_ref_local,
    N_ref_local,
    delta,
    M,
    time_steps,
    A,
    Ns,
    tau,
    thetas,
    coupled_A
  );
}


// [[Rcpp::export]]
Rcpp::List test_xpg_coupling_cpp(
    const arma::vec& Mu,
    double delta,
    int B
) {
  
  int M = Mu.n_elem - 1;
  
  arma::mat N_out(
    B,
    M + 1,
    arma::fill::zeros
  );
  
  arma::mat tilde_out(
    B,
    M + 1,
    arma::fill::zeros
  );
  
  arma::mat private_out(
    B,
    M + 1,
    arma::fill::zeros
  );
  
  
  // Sort means once because Mu is fixed
  arma::uvec sorted_col =
    arma::sort_index(Mu);
  
  Rcpp::IntegerVector sorted_indices(
    sorted_col.begin(),
    sorted_col.end()
  );
  
  
  // Position of reference particle 0
  int c =
    std::find(
      sorted_indices.begin(),
      sorted_indices.end(),
      0
    ) -
    sorted_indices.begin() + 1;
  
  
  for (int b = 0; b < B; ++b) {
    
    // ---------------------------------------------
    // Reference count
    // ---------------------------------------------
    
    double N_ref =
      R::rpois(Mu(0));
    
    N_out(b, 0) = N_ref;
    
    
    // ---------------------------------------------
    // Shared component of reference
    // ---------------------------------------------
    
    double tilde_ref =
      R::rbinom(
        N_ref,
        1.0 - delta
      );
    
    
    arma::vec tilde(
      M + 1,
      arma::fill::value(
        arma::datum::nan
      )
    );
    
    tilde(0) = tilde_ref;
    
    
    // ---------------------------------------------
    // Same coupling used in production
    // ---------------------------------------------
    
    if (c == 1) {
      
      tilde =
        update_proposals_cpp(
          make_sequence(c, M),
          "forward",
          sorted_indices,
          Mu,
          tilde,
          delta
        );
      
    } else if (c == M + 1) {
      
      tilde =
        update_proposals_cpp(
          make_sequence(c, 2),
          "backward",
          sorted_indices,
          Mu,
          tilde,
          delta
        );
      
    } else {
      
      tilde =
        update_proposals_cpp(
          make_sequence(c, 2),
          "backward",
          sorted_indices,
          Mu,
          tilde,
          delta
        );
      
      tilde =
        update_proposals_cpp(
          make_sequence(c, M),
          "forward",
          sorted_indices,
          Mu,
          tilde,
          delta
        );
    }
    
    
    tilde_out.row(b) =
      tilde.t();
    
    
    // ---------------------------------------------
    // Private components
    // ---------------------------------------------
    
    for (int j = 1; j < M + 1; ++j) {
      
      double private_j =
        R::rpois(
          delta * Mu(j)
        );
      
      private_out(b, j) =
        private_j;
      
      N_out(b, j) =
        tilde(j) +
        private_j;
    }
  }
  
  
  return Rcpp::List::create(
    Rcpp::Named("N") =
      N_out,
    
    Rcpp::Named("tilde") =
      tilde_out,
    
    Rcpp::Named("private") =
      private_out
  );
}


// [[Rcpp::export]]
Rcpp::List propagate_then_runif_cpp(
    int l,
    arma::mat X_ancestors,
    const arma::mat& X_ref_local,
    const arma::mat& N_ref_local,
    double delta,
    int M,
    const arma::vec& time_steps,
    const arma::mat& A,
    double N,
    double tau,
    const arma::vec& thetas
) {
  
  Rcpp::List prop =
    propagate_rcpp_SEIR_single(
      l,
      X_ancestors,
      X_ref_local,
      N_ref_local,
      delta,
      M,
      time_steps,
      A,
      N,
      tau,
      thetas
    );
  
  double U = R::runif(0.0, 1.0);
  
  return Rcpp::List::create(
    Rcpp::Named("propagation") = prop,
    Rcpp::Named("U") = U
  );
}


// [[Rcpp::export]]
Rcpp::List propagate_SIR_then_runif_cpp(
    int l,
    arma::mat X_ancestors,
    const arma::mat& X_ref_local,
    const arma::mat& N_ref_local,
    double delta,
    int M,
    const arma::vec& time_steps,
    const arma::mat& A,
    const arma::vec& Ns,
    double tau,
    const arma::vec& thetas,
    const arma::mat& coupled_A
) {

  Rcpp::List prop =
    propagate_rcpp_SIR(
      l,
      X_ancestors,
      X_ref_local,
      N_ref_local,
      delta,
      M,
      time_steps,
      A,
      Ns,
      tau,
      thetas,
      coupled_A
    );

  double U = R::runif(0.0, 1.0);

  return Rcpp::List::create(
    Rcpp::Named("propagation") = prop,
    Rcpp::Named("U") = U
  );
}





