//
// Created by Giulia de Sanctis on 01/12/24.
//

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    mat A(4, 5, fill::randu);
    mat B(4, 5, fill::randu);

    cout << A*B.t() << endl;

    return 0;
}

/*
int const q=5;
// dati se li carico direttamente da codice
//const int p= t[0].size()
//const int t=t.size()
//std::vector<double> tij, tg;
//    for (double t = 0; t <= 46; t += 2) {
//        tij.push_back(t / 46);
//    }
//    for (double t = 0; t <= 46; t += 0.1) {
//        tg.push_back(t / 46);
//    }
//    std::vector<std::vector<double>> B(T, std::vector<double>(2 * q, 0.0));
//    std::vector<std::vector<double>> Bpred(tg.size(), std::vector<double>(2 * q, 0.0));
//    std::vector<double> lambda2 = {8, 12, 16, 24, 48};
//    std::vector<double> periods(lambda2.size());
//    std::vector<double> lambda(lambda2.size());
//    for (size_t i = 0; i < lambda2.size(); ++i) {
//        periods[i] = lambda2[i] / 2;
//        lambda[i] = periods[i] / 46;
//    }


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

