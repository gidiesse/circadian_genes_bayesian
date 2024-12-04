//
// Created by Giulia de Sanctis on 01/12/24.
//

#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace std;
using namespace arma;

 void load_matrix (const std::string& file_path, const std::string& file_name, arma::mat &mat, bool print);
 void load_vector (const std::string& file_path, const std::string& file_name, arma::vec &vec, bool print); 
 void load_constant (const std::string& file_path, const std::string& file_name, int num, bool print);

int main()
{
    //arma::set_print_format(arma::PrintFormat::scientific); 
    std::string file_path = "/Users/giuliadesanctis/Desktop/POLIMI/mag_4_SEM_39/Bayesian/Progetto/Data/";
    
    // Define and load the matrices 
    arma::mat eta, theta, theta_tilde, thr1, W, Y; 
    // arma::mat B; 
    // arma::mat B_pred; 
    arma::mat lambda; 
    //load_matrix(file_path, "B.csv", B, true);
    //load_matrix(file_path, "Bpred.csv", Bpred, false);
    //load_matrix(file_path, "eta.csv", eta, false);
    //load_matrix(file_path, "Lambda.csv", lambda, false);
    //load_matrix(file_path, "THETA.csv", theta, false);
    //load_matrix(file_path, "Thetatilde.csv", theta_tilde, false);
    //load_matrix(file_path, "Thr1.csv", thr1, false);
    //load_matrix(file_path, "W.csv", W, false);
    //arma::mat Y; 
    load_matrix(file_path, "Y.csv", Y, false);

    // Define and load the vectors - no prob with numeff, probs with tg and tij for rounding
    //arma::vec numeff, tg, tij; 
    //load_vector(file_path, "numeff.csv", numeff, false);
    //load_vector(file_path, "tg.csv", tg, false);
    //load_vector(file_path, "tij.csv", tij, false);  

    // Define and load the constants
    //int p,T;
    //load_constant(file_path, "p.csv", p, false);
    //load_constant(file_path, "T.csv", T, false);

    // Case where I have real data and so only load Y 
    int p = Y.n_cols;
    int T = Y.n_rows;
    int const q=5;

    arma::vec t_ij = arma::linspace(0,46,24)/46;
    arma::vec tg = arma:: linspace(0, 46, 461)/46; 

    arma::mat B = zeros<mat>(T,2*q); 
    arma::mat B_pred = zeros<mat>(size(tg)[0], 2*q);
    arma::vec initial_periods = arma::vec("8, 12, 16, 24, 48");
    arma::vec periods = initial_periods / 2; 
    arma::vec full_periods = periods / 46; 

    for (size_t h=0; h<size(full_periods)[0]; h++) {
        arma::vec val1 = sin(((2*datum::pi)/(full_periods(h)))*t_ij);
        B.col(2*h) = val1; 
        arma::vec val2 = sin(((2*datum::pi)/(full_periods(h)))*tg);
        B_pred.col(2*h) = val2;
        arma::vec val3 = cos(((2*datum::pi)/(full_periods(h)))*t_ij);
        B.col(2*h+1) = val3;
        arma::vec val4 = cos(((2*datum::pi)/(full_periods(h)))*tg);
        B_pred.col(2*h+1) = val4;
    }

arma::arma_rng::set_seed(230518); // set the seed of the random number generator
// parameters gibbs sampler
int rep = 1; 
int nrun = 10000;         // Tot. number of iteration
int  burn = 1000;          // Burn-in
int thin = 5;             // Thinning
double sp = (nrun - burn)/thin; //Number of posterior samples
int k = 5;  //Number of factors to start with (for now)
// int k=int(log(p)*4); //Number of factors to start with (number of columns of Lambda)
                                   
double b0 = 1, b1 = 0.0005; //parameters for exponential - assumed to be arbitrary
double epsilon = 1e-3;  //Threshold limit
double prop = 1.00; //Proportion of redundant elements within columns

double as = 1, bs = 0.5; //Gamma hyperparameters for residual precision (true value res variance = 1 for every i)
double df = 3; //Gamma hyperparameters for t_{ij} [for phi_{ij} (i.e. rho)?]
double ad1 = 2.1, bd1 = 1; //Gamma hyperparameters for delta_1
double ad2 = 3.1, bd2 = 1; //gamma hyperparameters delta_h, h >= 2
double adf = 1, bdf = 1;          //Gamma hyperparameters for ad1 and ad2 or df

// Initial values
arma::colvec sig(p);

std::cout << "test 96"; 

for (int i=0; i<p; i++ ) {
    sig(i)=1/as_scalar(arma::randg(1, arma::distr_param(as,1/bs)));  //Residual variance (diagonal of sig^(-2)_i across proteins)
  }
 
arma::vec odd = arma::linspace(1,q,3);
arma::vec even = arma::linspace(0,q-1,3);

lambda = arma::zeros(p,k);  //loading matrix
eta =  arma::mvnrnd(arma::colvec(T,k, fill::zeros), arma::mat(k,k,fill::eye)); //Latent factors (distrib. = mean 0 and identity cov matrix)
std::cout << "\n AIUTOOOOOO";
W = arma::mvnrnd(arma::zeros(2*q,k), arma::eye(k,k)).t(); //Low dim. matrix W
theta_tilde = lambda * W;  //Matrix of unshrunk coefficients
theta = arma::zeros(p,2*q);  //Matrix of (fixed) basis functions coefficients
int Kappatheta = 5;
thr1 = randu( p, q, distr_param(0,Kappatheta)); //Matrix of thresholds for THETA



/*for (int i=0; i<q; i++) {
	arma::vec hypot_values = sqrt(pow(theta_tilde.col(2*i), 2) + pow(theta_tilde.col(2*i+1), 2));
    // Creare un vettore booleano per il confronto
    arma::uvec index = find(hypot_values >= thr1.col(i));
    // Se ci sono indici, aggiorna la matrice THETA
    if (index.n_elem > 0) {
        std::cout << index.n_elem; 
        theta(index, {2*i-1, 2*i}) = theta_tilde(index, {2*i-1, 2*i});
    }
    }*/ 
int i = 1;
std::cout << sqrt(pow(theta_tilde.col(2*i), 2) + pow(theta_tilde.col(2*i+1), 2));

    return 0;
}

void load_matrix (const std::string& file_path, const std::string& file_name, arma::mat &mat, bool print) {
        std::string specific_path = file_path + file_name;
        mat.load(specific_path, arma::csv_ascii);
        if (print == true)
            //std::cout << mat << endl; 
            mat.raw_print(); 
    }

void load_vector (const std::string& file_path, const std::string& file_name, arma::vec &vec, bool print) {
        std::string specific_path = file_path + file_name;
        vec.load(specific_path, arma::csv_ascii);
        if (print == true)
            vec.raw_print(); 
    }

void load_constant (const std::string& file_path, const std::string& file_name, int num, bool print) {
        std::string specific_path = file_path + file_name;
        std::ifstream num_file(specific_path);
         num_file >> num;
        if (print == true)
            std::cout << num << endl; 
    }


/*


for(size_t h=0;h<lambda.size();h++){
double lambda_h = lambda[h];
double factor = 2 * M_PI / lambda_h;
for (size_t i = 0; i < T; ++i) {
B[i][2 * h] = std::sin(factor * tij[i]);
Bpred[i][2 * h] = std::sin(factor * tg[i]);

B[i][2 * h + 1] = std::cos(factor * tij[i]);
Bpred[i][2 * h + 1] = std::cos(factor * tg[i]);
}
}
const size_t ktr = 20;
const size_t rep = 1;
const size_t nrun = 50;
const size_t burn = 20;
const size_t thin = 5;
const size_t sp = (nrun - burn)/thin; //Number of posterior samples
const size_t k = 5;  // Number of factors to start with (for now)

const size_t b0 = 1;
const size_t b1 = 0.0005;
const size_t epsilon = 1e-3; // Threshold limit
const size_t prop = 1.00; //Proportion of redundant elements within columns
// Define hyperparameter values
const size_t as = 1;
const size_t bs = 0.5;                           // Gamma hyperparameters for residual precision (true value res variance = 1 for every i)
const size_t df = 3;                                     // Gamma hyperparameters for t_{ij} ?? dubbio non sarebbe di phi
const size_t ad1 = 2.1;
const size_t bd1 = 1;                         // Gamma hyperparameters for delta_1
const size_t ad2 = 3.1;
const size_t bd2 = 1;                 // gamma hyperparameters delta_h, h >= 2
const size_t adf = 1;
const size_t bdf = 1;                           // Gamma hyperparameters for ad1 and ad2 or df

//initial values
std::gamma_distribution<double> gamma_dist(as, 1.0 / bs);
std::vector<double> sig(p);
for (int i = 0; i < p; ++i) {
sig[i] = gamma_dist(generator);
}
std::vector<int> odd;
std::vector<int> even;
for (size_t i=1;i<=q;i++){
odd.push_back(2*i-1);
even.push_back(2*i);
}
}
*/

