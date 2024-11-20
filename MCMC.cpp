#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>

using namespace std;
int main(){
    mt19937 generator(500);
    std::ifstream Yfile("Y.csv");
    std::string line;
    std::vector<std::vector<double>> Y;

    while (std::getline(Yfile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice

        // Leggi i valori separati da virgola
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        Y.push_back(row); // Aggiungi la riga alla matrice
    }
    Yfile.close();
    std::ifstream tijfile("tij.csv");
    std::vector<double> tij;

    while (std::getline(tijfile, line)) {
        std::stringstream ss(line);
        std::string value;
        tij.push_back(std::stod(value)); // Converte la stringa in double
    }
    tijfile.close();

    std::vector<double> tg;
    std::ifstream tgfile("tg.csv");
    while (std::getline(tgfile, line)) {
        std::stringstream ss(line);
        std::string value;
        tg.push_back(std::stod(value)); // Converte la stringa in double
    }
    tgfile.close();
    int p;
    std::ifstream pfile("p.csv");
    while (std::getline(pfile, line)) {
        std::stringstream ss(line);
        std::string value;
        p=std::stod(value); // Converte la stringa in double
    }
    pfile.close();
    int T;
    std::ifstream Tfile("T.csv");
    while (std::getline(Tfile, line)) {
        std::stringstream ss(line);
        std::string value;
        T=std::stod(value); // Converte la stringa in double
    }
    Tfile.close();
    std::vector<std::vector<double>> B;
    std::ifstream Bfile("B.csv");
    while (std::getline(Bfile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        B.push_back(row); // Aggiungi la riga alla matrice
    }
    Bfile.close();
    std::vector<std::vector<double>> eta;
    std::ifstream etafile("eta.csv");
    while (std::getline(etafile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        eta.push_back(row); // Aggiungi la riga alla matrice
    }
    etafile.close();
    std::vector<double> lambda;
    std::ifstream Lambdafile("Lambda.csv");

    if (std::getline(Lambdafile, line)) {
        std::stringstream ss(line); // Usa stringstream per analizzare la riga
        std::string value;
        while (std::getline(ss, value, ',')) {
            lambda.push_back(std::stoi(value)); // Converti la stringa in int e aggiungi al vettore
        }
    }
    Lambdafile.close();
    std::vector<std::vector<double>> thetatilde;
    std::ifstream Thetatildefile("Thetatilde.csv");
    while (std::getline(Thetatildefile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        thetatilde.push_back(row); // Aggiungi la riga alla matrice
    }
    Thetatildefile.close();
    std::vector<std::vector<double>> theta;
    std::ifstream THETAfile("THETA.csv");
    while (std::getline(THETAfile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        theta.push_back(row); // Aggiungi la riga alla matrice
    }
    THETAfile.close();
    std::vector<std::vector<double>> W;
    std::ifstream Wfile("W.csv");
    while (std::getline(Wfile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        W.push_back(row); // Aggiungi la riga alla matrice
    }
    Wfile.close();
    std::vector<double> numeff;
    std::ifstream numefffile("numeff.csv");
    if (std::getline(Lambdafile, line)) {
        std::stringstream ss(line); // Usa stringstream per analizzare la riga
        std::string value;
        while (std::getline(ss, value, ',')) {
            lambda.push_back(std::stoi(value)); // Converti la stringa in int e aggiungi al vettore
        }
    }
    numefffile.close();
    std::vector<std::vector<double>> Bpred;
    std::ifstream Bpredfile("Bpred.csv");
    while (std::getline(Bpredfile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        Bpred.push_back(row); // Aggiungi la riga alla matrice
    }
    Bpredfile.close();
    std::vector<std::vector<double>> thr1;
    std::ifstream thr1file("thr1.csv");
    while (std::getline(thr1file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        thr1.push_back(row); // Aggiungi la riga alla matrice
    }
    thr1file.close();
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
    const size_t sp = (nrun - burn)/thin;
    const size_t k = 5;

    const size_t b0 = 1;
    const size_t b1 = 0.0005;
    const size_t epsilon = 1e-3;
    const size_t prop = 1.00;
    // Define hyperparameter values
    const size_t as = 1;
    const size_t bs = 0.5;                           // Gamma hyperparameters for residual precision (true value res variance = 1 for every i)
    const size_t df = 3;                                     // Gamma hyperparameters for t_{ij}
    const size_t ad1 = 2.1;
    const size_t bd1 = 1;                         // Gamma hyperparameters for delta_1
    const size_t ad2 = 3.1;
    const size_t bd2 = 1;                 // gamma hyperparameters delta_h, h >= 2
    const size_t adf = 1;
    const size_t bdf = 1;                           // Gamma hyperparameters for ad1 and ad2 or df

//initial values
}


