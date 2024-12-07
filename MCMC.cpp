#include <iostream>
#include <armadillo>

void load_matrix (const std::string& file_path, const std::string& file_name, arma::mat &mat, bool print);
void load_vector (const std::string& file_path, const std::string& file_name, arma::vec &vec, bool print);
void load_constant (const std::string& file_path, const std::string& file_name, int num, bool print);

int main()
{
    std::string file_path = "D:/Politecnico/MAGISTRALE/Corsi/Bayesian statistics/Project/Check/Data/";
    arma::arma_rng::set_seed(230518);

    arma::mat Y;
    load_matrix(file_path, "Y.csv", Y, false);

    int p = Y.n_cols;
    int T = Y.n_rows;

    const int q = 5;

    arma::vec t_ij = arma::linspace<arma::vec>(0,46,24);
    t_ij = t_ij / arma::max(t_ij);    // standardized time points

    arma::vec tg = arma::linspace<arma::vec>(0, 46, 461) / 46;

    arma::mat B(T, 2*q, arma::fill::zeros);
    arma::mat B_pred(tg.n_elem, 2*q, arma::fill::zeros);

    arma::vec initial_periods = arma::vec("8, 12, 16, 24, 48");
    arma::vec periods = initial_periods / 2;
    arma::vec full_periods = periods / 46;

    for (size_t h = 0; h < full_periods.n_elem; ++h) {
        B.col(2*h) = sin(((2*arma::datum::pi) / (full_periods(h)))*t_ij);
        B_pred.col(2*h) = sin(((2*arma::datum::pi) / (full_periods(h)))*tg);

        B.col(2*h+1) = cos(((2*arma::datum::pi)/(full_periods(h)))*t_ij);
        B_pred.col(2*h+1) = cos(((2*arma::datum::pi)/(full_periods(h)))*tg);
    }

    // Define global constants
    int rep = 1;
    int nrun = 10000;                   // Tot. number of iteration
    int burn = 1000;                    // Burn-in
    int thin = 5;                       // Thinning
    double sp = (nrun - burn) / thin;   // Number of posterior samples
    int k = 5;                          // Number of factors to start with (for now)
    // int k = int(log(p)*4);           // Number of factors to start with (number of columns of Lambda)

    double b0 = 1, b1 = 0.0005;         // Parameters for exponential - assumed to be arbitrary
    double epsilon = 1e-3;              // Threshold limit
    double prop = 1.00;                 // Proportion of redundant elements within columns

    // Define hyper-parameter values
    double as = 1, bs = 0.5;            // Gamma hyper-parameters for residual precision (true value res variance = 1 for every i)
    double df = 3;                      // Gamma hyper-parameters for t_{ij} [for phi_{ij} (i.e. rho)?]
    double ad1 = 2.1, bd1 = 1;          // Gamma hyper-parameters for delta_1
    double ad2 = 3.1, bd2 = 1;          // Gamma hyper-parameters delta_h, h >= 2
    double adf = 1, bdf = 1;            // Gamma hyper-parameters for ad1 and ad2 or df

    // Initial values

    // scalar_type s = randg<scalar_type>( distr_param(a,b) ), where scalar_type is either float or double
    // a: shape parameter, b: scale parameter
    arma::colvec sig = arma::randg(p, arma::distr_param(as, 1/bs));

    arma::vec odd = arma::regspace(1, 2, q);
    arma::vec even = arma::regspace(0, 2, q-1);

    arma::mat lambda(p, k, arma::fill::zeros);      // Factor loadings
    arma::mat eta(T, k, arma::fill::zeros);         // Latent factors

    // TO BE CONTINUED

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
