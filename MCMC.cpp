//
// Created by Giulia de Sanctis on 01/12/24.
//

#include <iostream>
#include <armadillo>

void load_matrix (const std::string& file_path, const std::string& file_name, arma::mat &mat, bool print);
void load_vector (const std::string& file_path, const std::string& file_name, arma::vec &vec, bool print);
void load_constant (const std::string& file_path, const std::string& file_name, int num, bool print);

int main()
{
    std::string file_path = "../Data/";
    arma::arma_rng::set_seed(230518);

    arma::mat Y;
    load_matrix(file_path, "Y.csv", Y, false);

    int p = Y.n_cols;       // p = number of proteins
    int T = Y.n_rows;       // T = number of time points

    const int q = 5;        // number of sin (cos) bases, for a total of 2q bases - no intercept

    arma::vec t_ij = arma::linspace<arma::vec>(0, 46, 24);
    t_ij = t_ij / arma::max(t_ij);                              // standardized time points

    arma::vec tg = arma::linspace<arma::vec>(0, 46, 461) / 46;  // standardized prediction grid

    arma::mat B(T, 2*q, arma::fill::zeros);                     // design matrix (data)
    arma::mat B_pred(tg.n_elem, 2*q, arma::fill::zeros);        // design matrix (estimation/prediction)

    arma::rowvec initial_periods = {8, 12, 16, 24, 48};         // range of periods for fixed basis functions
    arma::rowvec periods = initial_periods / 2;                 // periods on a unit-based increment
    arma::rowvec full_periods = periods / 46;                   // periods for our time increments

    for (size_t h = 0; h < full_periods.n_elem; ++h)
    {
        B.col(2*h) = sin(((2*arma::datum::pi) / (full_periods(h)))*t_ij);
        B_pred.col(2*h) = sin(((2*arma::datum::pi) / (full_periods(h)))*tg);

        B.col(2*h+1) = cos(((2*arma::datum::pi)/(full_periods(h)))*t_ij);
        B_pred.col(2*h+1) = cos(((2*arma::datum::pi)/(full_periods(h)))*tg);
    }

    // Define global constants
    int rep = 1;
    int nrun = 10000;                   // Number of iteration
    int burn = 1000;                    // Burn-in
    int thin = 5;                       // Thinning
    double sp = (nrun - burn) / thin;   // Number of posterior samples
    int k = 5;                          // Number of factors to start with (for now)
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

    // Initial values

    // scalar_type s = randg<scalar_type>( distr_param(a,b) ), where scalar_type is either float or double
    // a: shape parameter, b: scale parameter
    arma::colvec sig = arma::randg(p, arma::distr_param(as, 1/bs));

    arma::rowvec odd = arma::linspace<arma::rowvec>(1, 2*q - 1, q);
    arma::rowvec even = arma::linspace<arma::rowvec>(2, 2*q, q);

    arma::mat lambda(p, k, arma::fill::zeros);              // Factor loadings

    arma::mat eta = arma::mvnrnd(arma::colvec(k, arma::fill::zeros),
                                 arma::eye(k, k), T).t();   // Latent factors

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

    arma::colvec tau_h = arma::cumprod(delta);                             // Global shrinkage coefficients

    arma::mat P_lam(p,k);                            // Precision of loadings rows
    for (size_t i = 0; i < p; ++i)
        P_lam.row(i) = phi_ih.row(i) % tau_h.t();

    arma::mat acc2(nrun, p, arma::fill::zeros);      // Matrix of rejection probabilities of MH steps


    // -----  Gibbs sampler  -----

    for (size_t i = 0; i < 1; ++i) {

        // Update error precisions
        arma::mat Y_til = Y - B*theta.t() - eta*lambda.t();

        double a = as + 0.5 * T;
        arma::mat b = arma::ones(1, p) / (bs + 0.5 * arma::sum(arma::pow(Y_til, 2)));

        arma::vec sig(p);
        for (size_t j = 0; j < p; ++j)
            sig(j) = arma::randg(arma::distr_param(a,b(j)));


        // Update eta
        arma::mat lmsg(p,k);
        for (size_t j = 0; j < k; ++j)
            lmsg.col(j) = lambda.col(j) % sig;
        arma::mat veta_1 = arma::mat(k,k, arma::fill::eye) + lmsg.t() * lambda;

        // TO BE CONTINUED

        
    }

    return 0;
}

void load_matrix (const std::string& file_path, const std::string& file_name, arma::mat &mat, bool print=true) {
    std::string specific_path = file_path + file_name;
    mat.load(specific_path, arma::csv_ascii);
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
