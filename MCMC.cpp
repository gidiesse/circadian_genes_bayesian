#include <iostream>
#include <armadillo>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <RcppArmadillo.h>

void load_matrix (const std::string& file_path, const std::string& file_name, arma::mat &mat, bool print);
void load_vector (const std::string& file_path, const std::string& file_name, arma::vec &vec, bool print);
void load_constant (const std::string& file_path, const std::string& file_name, int num, bool print);
arma::uvec setdiff(const arma::vec& A, const arma::vec& B);
void write_matrix (const std::string& file_path, arma::mat &mat,const std::string& name_csv);
arma::mat standardizeMatrix(const arma::mat& matrix); 
arma::mat logTransformation(const arma::mat& inputMatrix);


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void MCMC_Circadian_Genes(arma::mat Yn, arma::vec t_ij, arma::vec tg, const std::string& file_path,const int nrun,std::vector<std::string> vars_to_save,int burn)
{
  auto start=std::chrono::high_resolution_clock::now();
  arma::arma_rng::set_seed(250);
  Yn=Yn.t();
  arma::mat Y1;
  arma::mat Y;

  Y=standardizeMatrix(Yn);

  int max_latent_factors = 0;
  
  int T = Y.n_rows;       // T = number of time points
  int p = Y.n_cols;       // p = number of proteins
  
  const int q = 5;        // number of sin (cos) bases, for a total of 2q bases - no intercept

  t_ij = t_ij / arma::max(t_ij);                              // standardized time points

  arma::mat B(T, 2*q, arma::fill::zeros);                     // design matrix (data)
  arma::mat B_pred(tg.n_elem, 2*q, arma::fill::zeros);        // design matrix (estimation/prediction)
  
  arma::rowvec initial_periods = {8, 12, 16, 20, 24};         // range of periods for fixed basis functions
  arma::rowvec full_periods = initial_periods / 48;
  for (size_t h = 0; h < full_periods.n_elem; ++h)
  {
        B.col(2*h) = sin(((2*arma::datum::pi) / (full_periods(h)))*t_ij);
        B_pred.col(2*h) = sin(((2*arma::datum::pi) / (full_periods(h)))*tg);
        B.col(2*h+1) = cos(((2*arma::datum::pi)/(full_periods(h)))*t_ij);
        B_pred.col(2*h+1) = cos(((2*arma::datum::pi)/(full_periods(h)))*tg);
  }
  // Define global constants
  int rep = 1;
  //int nrun = 1000;                   // Number of iteration
  // int burn = 100;                    // Burn-in
  int thin = 5;                       // Thinning
  double sp = (nrun - burn) / thin;   // Number of posterior samples
  double k = 5;                          // Number of factors to start with (for now)
  // int k = arma::floor(log(p)*4);   // Number of factors to start with (number of columns of Lambda)

  double b0 = 1., b1 = 0.0005;        // Parameters for exponential - assumed to be arbitrary
  double epsilon = 1e-3;              // Threshold limit
  double prop = 1.00;                 // Proportion of redundant elements within columns

  // Define hyper-parameter values
  double as = 1., bs = 0.5;            // Gamma hyper-parameters for residual precision (true value res variance = 1 for every i)
  double df = 3.;                      // Gamma hyper-parameters for t_{ij} [for phi_{ij} (i.e. rho)?]
  double ad1 = 2.1, bd1 = 1.;          // Gamma hyper-parameters for delta_1
  double ad2 = 3.1, bd2 = 1.;          // Gamma hyper-parameters delta_h, h >= 2
  double adf = 1., bdf = 1.;           // Gamma hyper-parameters for ad1 and ad2 or df

  arma::mat Lambdaout(p * 12, sp, arma::fill::zeros); // 12 is the max number of latent factors
  arma::mat Etaout(T * 12, sp, arma::fill::zeros);
  arma::rowvec Thetaout(2 * q * p, arma::fill::zeros);
  arma::mat Thetaout_tot; //(sp, 2*q*p, arma::fill::zeros);

  // Initial values

  // scalar_type s = randg<scalar_type>( distr_param(a,b) ), where scalar_type is either float or double
  // a: shape parameter, b: scale parameter
  arma::colvec sig = arma::randg(p, arma::distr_param(as, 1/bs));

  arma::uvec odd = arma::linspace<arma::uvec>(1, 2*q - 1, q);
  arma::uvec even = arma::linspace<arma::uvec>(0, 2*q - 2, q);

  arma::mat lambda(p, k, arma::fill::zeros);              // Factor loadings

  arma::mat eta = arma::mvnrnd(arma::colvec(k, arma::fill::zeros),arma::eye(k, k), T).t();   // Latent factors

  arma::mat W = arma::mvnrnd(arma::colvec(k, arma::fill::zeros),
                               arma::eye(k, k), 2*q);

  arma::mat theta_tilde = lambda * W;                 // Matrix of un-shrunk coefficients
  arma::mat theta(p, 2*q, arma::fill::zeros);         // Matrix of (fixed) basis functions coefficients

  double kappa_theta = 5.;                            // Upper bound on the Uniform prior on the thresholds
  arma::mat thresholds = arma::randu(p, q, arma::distr_param(0., kappa_theta));  // Matrix of thresholds for theta

  for (size_t i = 0; i < q; ++i)
  {
        arma::uvec index = arma::find(arma::norm(theta_tilde.cols(2*i, 2*i+1), 2) >= thresholds.col(i));
        for(const auto& ind : index)
        {
            theta(ind, 2*i) = theta_tilde(ind, 2*i);
            theta(ind, 2*i + 1) = theta_tilde(ind, 2*i + 1);
        }
  }

  arma::mat phi_ih = arma::randg(p, k, arma::distr_param(df/2, 2/df));   // Local shrinkage coefficients

  arma::colvec delta(k, arma::fill::zeros);                              // Global shrinkage coefficients multipliers
  delta.row(0) = arma::randg(arma::distr_param(ad1,bd1));
  delta.rows(1,k-1) = arma::randg(k-1, 1, arma::distr_param(ad2,bd2));

  arma::colvec tau_h = arma::cumprod(delta);                      // Global shrinkage coefficients

  arma::mat P_lam(p,k);                            // Precision of loadings rows
  for (size_t i = 0; i < p; ++i)
        P_lam.row(i) = phi_ih.row(i) % tau_h.t();

  arma::mat acc2(nrun, p, arma::fill::zeros);      // Matrix of rejection probabilities of MH steps


    // -----  Gibbs sampler  -----

  for (size_t i = 0; i < nrun; ++i) {
    auto start_it=std::chrono::high_resolution_clock::now();

        // Update error precisions
    arma::mat Y_til = Y - B*theta.t() - eta*lambda.t();

    double a = as + 0.5 * T;
    arma::mat b = arma::ones(1, p) / (bs + 0.5 * arma::sum(arma::pow(Y_til, 2)));
    
    for (size_t j = 0; j < p; ++j) {
      sig(j) = arma::randg(arma::distr_param(a, b(j)));
    }
    
    // Update eta
    arma::mat lmsg(p,k);
    for (size_t j = 0; j < k; ++j) {
      lmsg.col(j) = lambda.col(j) % sig;
    }
    arma::mat V_eta_1 = arma::eye(k,k) + lmsg.t() * lambda;
    arma::mat T_chol_eta = arma::chol(V_eta_1);
    arma::mat Q_eta, R_eta;
    arma::qr(Q_eta, R_eta, T_chol_eta);
    arma::mat S_eta = arma::inv(R_eta);
    arma::mat V_eta = S_eta * S_eta.t();
    arma::mat M_eta = (( Y - B * theta.t()) * lmsg) * V_eta;
    eta = M_eta + arma::randn(T, k, arma::distr_param(0.,1.)) * S_eta.t();
    
    
    // Update of Lambda (Rue and Held - 2005)
    for (size_t h = 0; h < p; ++h) {
      arma::mat V_lambda_1 = sig(h) * eta.t() * eta + W * arma::eye(2*q, 2*q) * W.t() + arma::diagmat(P_lam.row(h));
      arma::mat T_chol_lambda = arma::chol(V_lambda_1);
      arma::mat Q_lambda, R_lambda;
      arma::qr(Q_lambda, R_lambda, T_chol_lambda);
      arma::mat S_lambda = arma::inv(R_lambda);
      arma::mat V_lambda= S_lambda * S_lambda.t();
      arma::rowvec M_lambda = (sig(h) * eta.t() * (Y.col(h) - B * theta.row(h).t()) +
        W * arma::eye(2*q, 2*q) * theta_tilde.row(h).t()).t() * V_lambda;
      lambda.row(h) = M_lambda + arma::randn(1, k, arma::distr_param(0.,1.)) * S_lambda.t();
    }
    
    // Update phi_ih
    arma::mat b_den_prod(p,k);
    for (size_t j = 0; j < p; ++j) {
      b_den_prod.row(j) = arma::pow(lambda.row(j),2) % tau_h.t();
    }
    a = df/2 + 0.5;
    b = arma::ones(p, k) / (df/2 + b_den_prod);
    for (size_t j = 0; j < p; ++j) {
      for (size_t l = 0; l < k; ++l) {
        phi_ih(j,l) = arma::randg(arma::distr_param(a,b(j,l)));
      }
    }
    
    
    // Update delta and tau_h
    arma::mat b_mat(p,k);
    for (size_t j = 0; j < p; ++j) {
      b_mat.row(j) = arma::pow(lambda.row(j),2) % phi_ih.row(j);
    }
    double ad = ad1 + 0.5 * p * k;
    double bd = bd1 + 0.5 * (1/delta(0)) * arma::sum(tau_h % arma::sum(b_mat).t());
    delta(0) = arma::randg(arma::distr_param(ad,1/bd));
    tau_h = arma::cumprod(delta);
    
    for (size_t h = 1; h < k; ++h) {
      ad = ad2 + 0.5 * p * (k-h);
      bd = bd2 + 0.5 * (1/delta(h)) * arma::sum(tau_h(arma::span(h,tau_h.n_elem-1)) % arma::sum(b_mat.cols(h,arma::size(b_mat)[1]-1)).t());
      delta(h) = arma::randg(arma::distr_param(ad, 1/bd));
      tau_h = arma::cumprod(delta);
    }
    
    
    // Update precision parameters
    for (size_t j = 0; j < p; ++j) {
      P_lam.row(j) = phi_ih.row(j) % tau_h.t();
    }
    
    // Update matrix W
    arma::mat V_W_1 = lambda.t() * lambda + arma::eye(k,k);
    arma::mat T_chol_W = arma::chol(V_W_1);
    arma::mat Q_W, R_W;
    arma::qr(Q_W, R_W, T_chol_W);
    arma::mat S_W = arma::inv(R_W);
    arma::mat V_W = S_W * S_W.t();
    
    for (size_t h = 0; h < 2*q; ++h) {
      arma::mat M_W = (lambda.t() * theta_tilde.col(h)).t() * V_W;
      W.col(h) = (M_W + arma::randn(1,k, arma::distr_param(0.,1.)) * S_W.t()).t();
    }
    
    
    // Update theta_tilde
    for (size_t h = 0; h < p; ++h) {
      arma::mat V_theta_tilde_1 = sig(h) * B.t() * B + arma::eye(2*q,2*q);
      arma::mat T_chol_theta_tilde = arma::chol(V_theta_tilde_1);
      arma::mat Q_theta_tilde, R_theta_tilde;
      arma::qr(Q_theta_tilde, R_theta_tilde, T_chol_theta_tilde);
      arma::mat S_theta_tilde = arma::inv(R_theta_tilde);
      arma::mat V_theta_tilde = S_theta_tilde * S_theta_tilde.t();
      
      arma::vec Y_cent = Y.col(h) - eta * lambda.row(h).t();
      arma::mat W_trans = W.t();
      arma::rowvec M_theta_tilde = (sig(h) * B.t() * Y_cent + W_trans * lambda.row(h).t()).t() * V_theta_tilde;
      arma::rowvec theta_tilde_prop = M_theta_tilde + arma::randn(1,2*q,arma::distr_param(0.,1.)) * S_theta_tilde.t();
      arma::rowvec theta_prop(2*q, arma::fill::zeros);
      
      arma::uvec index = arma::find(arma::sqrt(arma::pow(theta_tilde_prop.elem(odd).t(),2) + arma::pow(theta_tilde_prop.elem(even).t(),2)) >= thresholds.row(h));
      
      if (index.n_elem != 0) {
        arma::uvec odd_index = odd(index);
        arma::uvec even_index = even(index);
        arma::uvec ind = arma::sort(arma::join_cols(odd_index, even_index));
        theta_prop(ind) = theta_tilde_prop(ind);
      }
      
      double r = arma::as_scalar(arma::exp( -0.5 * (theta_tilde_prop.t() - W_trans * lambda.row(h).t()).t() *
                                 (theta_tilde_prop.t() - W_trans * lambda.row(h).t()) +
                                 0.5 * (theta_tilde.row(h).t() - W_trans*lambda.row(h).t()).t() *
                                 (theta_tilde.row(h).t() - W_trans*lambda.row(h).t()) -
                                 0.5 * sig(h) * (Y.col(h) - B * theta_prop.t() - eta * lambda.row(h).t()).t() *
                                 (Y.col(h) - B * theta_prop.t() - eta * lambda.row(h).t()) +
                                 0.5 * sig(h) * (Y.col(h) - B * theta.row(h).t() - eta * lambda.row(h).t()).t() *
                                 (Y.col(h) - B * theta.row(h).t() - eta * lambda.row(h).t()) -
                                 0.5 * (theta_tilde.row(h).t() - M_theta_tilde.t()).t() *
                                 V_theta_tilde_1 * (theta_tilde.row(h).t() - M_theta_tilde.t()) +
                                 0.5 * (theta_tilde_prop.t() - M_theta_tilde.t()).t() *
                                 V_theta_tilde_1 * (theta_tilde_prop.t() - M_theta_tilde.t())));
      
      bool u = arma::as_scalar(arma::randu(1)) > std::fmin(1,r);
      acc2(i,h) = u;
      
      arma::rowvec theta_acc = theta_tilde_prop - (theta_tilde_prop - theta_tilde.row(h)) * u;
      theta_tilde.row(h) = theta_acc;
      theta.row(h) = arma::zeros(1, 2*q);
      
      arma::rowvec row_h = theta_tilde.row(h);
      index = arma::find(arma::sqrt(arma::pow(row_h.elem(odd).t(),2) + arma::pow(row_h.elem(even).t(),2)) >= thresholds.row(h));
      if (index.n_elem != 0) {
        arma::uvec odd_index = odd(index);
        arma::uvec even_index = even(index);
        arma::uvec ind = arma::sort(arma::join_cols(odd_index, even_index));
        for(const auto& id : ind) {
          theta(h,id) = theta_tilde(h,id);
        }
      }
    }
    
    
    // Update of thresholds of Theta
    for (size_t h = 0; h < q; ++h) {
      arma::vec v1 = arma::regspace(0,1,2*q-1);
      arma::vec v2 = arma::regspace(2*h,1,2*h+1);
      
      arma::vec AA = arma::vectorise(B.cols(setdiff(v1,v2)) * theta.cols(setdiff(v1,v2)).t());
      arma::vec BB = AA + arma::vectorise(B.cols(2*h,2*h+1) * theta_tilde.cols(2*h,2*h+1).t());
      
      arma::vec comp1 = arma::sqrt(arma::pow(theta_tilde.col(2*h),2) + arma::pow(theta_tilde.col(2*h+1),2));
      arma::vec Yprime= arma::vectorise(Y - eta*lambda.t());
      arma::vec nn = arma::diagvec(arma::reshape(Yprime - BB,T,p).t() * arma::reshape(Yprime - BB,T,p));
      arma::vec num = arma::exp(-0.5 * (sig % nn)) % comp1;
      arma::vec dd = arma::diagvec(arma::reshape(Yprime - AA,T,p).t() * arma::reshape(Yprime - AA,T,p));
      arma::vec den = num + arma::exp(-0.5 * (sig % dd)) % (kappa_theta - comp1);
      arma::rowvec short2 = (num / den).t();
      arma::uvec pc = arma::find((comp1 > kappa_theta));
      if (pc.n_elem != 0) {
        for (const auto& pp : pc) {
          thresholds(pp,h) = arma::randu(arma::distr_param(0.,kappa_theta));
        }
      }
      
      v1 = arma::regspace(0,1,p-1);
      arma::uvec npc = setdiff(v1,arma::conv_to<arma::vec>::from(pc));
      
      arma::rowvec u = arma::randu(1, short2.n_elem);
      for (const auto& pp : pc) {
        u(pp) = 0;
      }
      
      for (const auto& pp : npc) {
        // differentiate case where comp1(pp) = 0 since Armadillo syntax does not allow generation from uniform (0,0)
        if (comp1(pp) != 0)
          thresholds(pp,h) = (u(pp) <= short2(pp)) * arma::randu(arma::distr_param(0.,comp1(pp))) + (u(pp) > short2(pp)) * arma::randu(arma::distr_param(comp1(pp),kappa_theta));
        else
          thresholds(pp,h) = (u(pp) > short2(pp)) * arma::randu(arma::distr_param(comp1(pp),kappa_theta));
      }
      
      theta.cols(2*h,2*h+1) = arma::zeros(p,2);
      arma::uvec index = arma::find(comp1 >= thresholds.col(h));
      
      if (index.n_elem != 0) {
        for (const auto& pp : index) {
          //
          theta(pp, arma::span(2*h,2*h+1)) = theta_tilde(pp, arma::span(2*h, 2*h+1));
        }
      }
    }
    
    
    // Make adaptations
    // This part adapts the number of latent factors at the end of each run of the MCMC
    
    double prob = 1 / std::exp(b0 + b1 * i);
    double uu = arma::as_scalar(arma::randu(1));
    arma::rowvec lind = arma::conv_to<arma::rowvec>::from(arma::sum(arma::conv_to<arma::Mat<double>>::from(arma::abs(lambda) < epsilon), 0)) / static_cast<double>(p);
    arma::urowvec vettor = (lind >= prop);
    double num = arma::sum(vettor);
    
    if (uu < prob) {
      if (i > 20 && num == 0 && arma::all(lind < 0.995)) {
        ++k;
        lambda.insert_cols(k-1, arma::zeros(p,1));
        eta.insert_cols(k-1, arma::randn(T,1));
        W.insert_rows(k-1, arma::randn(1,2*q));
        phi_ih.insert_cols(k-1, arma::randg(p, 1, arma::distr_param(df/2, 2/df)));
        delta.resize(k);
        delta(k-1) = arma::randg(arma::distr_param(ad2,1/bd2));
        tau_h = arma::cumprod(delta);
        P_lam.insert_cols(k-1, arma::zeros(p,1));
        for (size_t j = 0; j < p; ++j) {
          P_lam.row(j) = phi_ih.row(j) % tau_h.t();
        }
        
      }
      else if (num > 0) {
        arma::vec v1 = arma::regspace(0,1,k-1);
        arma::vec v2 = arma::conv_to<arma::vec>::from(arma::find(vettor));
        arma::uvec non_red = setdiff(v1,v2);
        k = std::max(static_cast<double>(k)-num, 1.);
        lambda = lambda.cols(non_red);
        W = W.rows(non_red);
        phi_ih = phi_ih.cols(non_red);
        eta = eta.cols(non_red);
        delta = delta(non_red);
        tau_h = arma::cumprod(delta);
        P_lam = arma::zeros(p,k);
        for (size_t j = 0; j < p; ++j) {
          P_lam.row(j) = phi_ih.row(j) % tau_h.t();
        }
      }
    }
    
    if (k > max_latent_factors)
      max_latent_factors = k;
    
    std::cout << "number of latent factors is: " << k << std::endl;
    auto end_it=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> it_time = end_it - start_it;
    
    // Save sampled values (after thinning)
    if (i % thin == 0 && i>burn ) {
      // Calcola l'indice della colonna da aggiornare
      int col_index = (i-burn) / thin;
      
      arma::colvec add_col = arma::join_vert(arma::vectorise(lambda), arma::zeros<arma::colvec>(p * 12 - p * k));
      Lambdaout.col(col_index) = add_col;
      
      arma::colvec add_col2 = arma::join_vert(arma::vectorise(eta), arma::zeros<arma::colvec>(T * 12 - T * k));
      Etaout.col(col_index) = add_col2;
      
      Thetaout = arma::reshape(theta,1,p*2*q);
      Thetaout_tot = arma::join_vert(Thetaout_tot, Thetaout);
    }
    
    
  }
  
  
  std::cout << "Max number of latent factors reached: " << max_latent_factors << std::endl;
  
  std::unordered_map<std::string, arma::mat> map_matrix;
  //std::unordered_map<std::string, arma::vec> map_vec;
  map_matrix["theta"] = theta;
  map_matrix["theta_tilde"] = theta_tilde;
  map_matrix["eta"] = eta;
  map_matrix["W"] = W;
  map_matrix["B_pred"] = B_pred;
  map_matrix["B"] = B;
  map_matrix["Lambdaout"] = Lambdaout;
  map_matrix["Thetaout"] = Thetaout;
  map_matrix["Thetaout_tot"] = Thetaout_tot;
  map_matrix["Etaout"] = Etaout;
  map_matrix["thresholds"]=thresholds;
  
  
  for (std::size_t i = 0; i < vars_to_save.size(); ++i) {
    auto it = map_matrix.find(vars_to_save[i]);
    if (it == map_matrix.end()) {
      std::cerr << "The variable you are looking for doesn't exist" << std::endl;
    }
    else
      write_matrix(file_path,it->second,it->first);
  }
  //write_matrix (file_path,W);
  //write_matrix (file_path,lambda);
  //write_matrix (file_path,eta);
  //write_matrix (file_path,tau_h);
  //write_matrix (file_path,B);
  //write_matrix (file_path,B_pred);
  
  auto final=std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff=final-start;
  std::cout << diff.count() << std::endl;
  std::cout << "Max number of latent factors reached: " << max_latent_factors << std::endl;
}

// [[Rcpp::export]]
arma::uvec setdiff(const arma::vec& A, const arma::vec& B) {
  // Ottieni solo gli elementi unici e ordinati di A e B
  arma::vec uniqueA = arma::unique(arma::sort(A));
  arma::vec uniqueB = arma::unique(arma::sort(B));
  
  // Vettore per contenere gli elementi di A che non sono in B
  arma::uvec diff;
  
  // Trova gli elementi di uniqueA che non sono in uniqueB
  for (const auto& elem : uniqueA) {
    if (!arma::any(uniqueB == elem)) {
      diff.insert_rows(diff.n_rows, arma::uvec({static_cast<unsigned long long> (elem)}));
    }
  }
  
  return diff;
}

void load_matrix (const std::string& file_path, const std::string& file_name, arma::mat &mat, bool print=true) {
  std::string specific_path = file_path + file_name;
  mat.load(specific_path, arma::csv_ascii);
  mat.shed_row(0);    // Rimuove la prima riga (indice 0)
  if (print == true)
    mat.print();
}

void load_vector (const std::string& file_path, const std::string& file_name, arma::vec &vec, bool print=true) {
  std::string specific_path = file_path + file_name;
  vec.load(specific_path, arma::csv_ascii);
  
  if (print == true)
    vec.raw_print();
}

void load_constant (const std::string& file_path, const std::string& file_name, int num, bool print=true) {
  std::string specific_path = file_path + file_name;
  std::ifstream num_file(specific_path);
  num_file >> num;
  if (print == true)
    std::cout << num << std::endl;
}
// [[Rcpp::export]]
arma::mat standardizeMatrix(const arma::mat& matrix) {
  arma::mat standardized = matrix;
  
  for (size_t col = 0; col < matrix.n_cols; ++col) {
    double mean_value = arma::mean(matrix.col(col));
    double stddev_value = arma::stddev(matrix.col(col));
    
    if (stddev_value != 0) {
      standardized.col(col) = (matrix.col(col) - mean_value) / stddev_value;
    } else {
      std::cout << "Attention standard deviation is null in the column " << col <<std::endl;
    }
  }
  
  return standardized;
}
// [[Rcpp::export]]
void write_matrix (const std::string& file_path, arma::mat& mat, const std::string& name_csv) {
#ifdef _WIN32
  std::string separator = "\\";
#else
  std::string separator = "/";
#endif
  std::string full_path = file_path + separator + name_csv + ".csv";
  std::cout<<full_path<<std::endl;
  std::ofstream to_write(full_path);
  
  for (size_t i = 0; i < mat.n_rows; ++i) {
    for (size_t j = 0; j < mat.n_cols; ++j) {
      to_write << mat(i,j);
      if (j < mat.n_cols - 1) {
        to_write << ",";
      }
    }
    to_write << std::endl;
  }
  to_write.close();
  return;
}

// [[Rcpp::export]]
arma::mat logTransformation(const arma::mat& inputMatrix) {
  arma::mat outputMatrix = inputMatrix;  // Copia della matrice originale
  for (arma::uword i = 0; i < inputMatrix.n_rows; ++i) {
    for (arma::uword j = 0; j < inputMatrix.n_cols; ++j) {
      if (inputMatrix(i, j) > 0) {
        outputMatrix(i, j) = std::log(inputMatrix(i, j));  // Applicazione logaritmo naturale
      } else {
        std::cerr << "Errore: elemento non positivo in posizione (" << i << ", " << j << ").\n";
        outputMatrix(i, j) = arma::datum::nan;  // Assegna NaN se il valore è non positivo
      }
    }
  }
  return outputMatrix;
}
