#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>

using namespace std;

std::vector<std::vector<double>> read_matrix (std::ifstream file, std::vector<std::vector<double>> mat);
std::vector<double> read_vector (const std::string& file_path, std::vector<double> &vec);
void read_synthetic_data (std::string file_path);

int main() {
    mt19937 generator(500);
    std::string file_path = "/Users/giuliadesanctis/circadian_genes/circadian_genes_bayesian/Data/";
    read_synthetic_data(file_path);
}

std::vector<std::vector<double>> read_matrix (const std::string& file_path, std::vector<std::vector<double>> &mat) {
    std::string line;
    std::ifstream file(file_path);
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row; // Riga della matrice
// Leggi i valori separati da virgola
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value)); // Converte la stringa in double
        }
        mat.push_back(row); // Aggiungi la riga alla matrice
    }
    return mat;
}

std::vector<double> read_vector (const std::string& file_path, std::vector<double> &vec) {
    std::string line;
    std::ifstream file(file_path);
// Leggi i valori separati da virgola
        while (std::getline(file, line)) {
            try {
                vec.push_back(std::stod(line)); // Convert line to double and add to vector
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid line: '" << line << "' (not a valid number)\n";
            }
        }
    return vec;
}

void read_synthetic_data (std::string base_path) {
    typedef std::vector<std::vector<double>> matrix;
    typedef std::vector<double> vector;

    // read matrix Y
    matrix Y;
    std::string file_name = "Y.csv";
    std::string full_path = base_path + file_name;
    std::ifstream Y_file(full_path);
    // Check that it has opened correctly
    if (!Y_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open Y.csv\n";
    }
    Y = read_matrix(full_path, Y);
    for (auto i = 0; i <= 10; i ++) {
        for (auto j = 0; j <= 10; j++) {
            std::cout << Y[i][j] << "\t";
        }
        std::cout << "\n";
    }
    Y_file.close();

    // read matrix B
    matrix B;
    file_name = "B.csv";
    full_path = base_path + file_name;
    std::ifstream B_file(full_path);
    // Check that it has opened correctly
    if (!B_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open B.csv\n";
    }
    B = read_matrix(full_path, B);
    std::cout << "We print matrix B \n";
    for (auto i = 0; i <= 10; i ++) {
        for (auto j = 0; j <= 10; j++) {
            std::cout << B[i][j] << "\t";
        }
        std::cout << "\n";
    }
    B_file.close();

    // read matrix eta
    matrix eta;
    file_name = "eta.csv";
    full_path = base_path + file_name;
    std::ifstream eta_file(full_path);
    // Check that it has opened correctly
    if (!eta_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open eta.csv\n";
    }
    eta = read_matrix(full_path, eta);
    std::cout << "We print matrix eta \n";
    for (auto i = 0; i <= 5; i ++) {
        for (auto j = 0; j <= 5; j++) {
            std::cout << eta[i][j] << "\t";
        }
        std::cout << "\n";
    }
    eta_file.close();

    // read matrix Lambda
    matrix lambda;
    file_name = "Lambda.csv";
    full_path = base_path + file_name;
    std::ifstream lambda_file(full_path);
    // Check that it has opened correctly
    if (!lambda_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open lambda.csv\n";
    }
    lambda = read_matrix(full_path, lambda);
    std::cout << "We print matrix lambda \n";
    for (auto i = 0; i <= 5; i ++) {
        for (auto j = 0; j <= 5; j++) {
            std::cout << lambda[i][j] << "\t";
        }
        std::cout << "\n";
    }
    lambda_file.close();

    // read matrix theta_tilde
    matrix theta_tilde;
    file_name = "Thetatilde.csv";
    full_path = base_path + file_name;
    std::ifstream theta_tilde_file(full_path);
    // Check that it has opened correctly
    if (!theta_tilde_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open Thetatilde.csv\n";
    }
    theta_tilde = read_matrix(full_path, theta_tilde);
    std::cout << "We print matrix theta_tilde \n";
    for (auto i = 0; i <= 10; i ++) {
        for (auto j = 0; j <= 10; j++) {
            std::cout << theta_tilde[i][j] << "\t";
        }
        std::cout << "\n";
    }
    theta_tilde_file.close();

    // read matrix theta
    matrix theta;
    file_name = "Theta.csv";
    full_path = base_path + file_name;
    std::ifstream theta_file(full_path);
    // Check that it has opened correctly
    if (!theta_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open THETA.csv\n";
    }
    theta = read_matrix(full_path, theta);
    std::cout << "We print matrix theta \n";
    for (auto i = 0; i <= 10; i ++) {
        for (auto j = 0; j <= 10; j++) {
            std::cout << theta[i][j] << "\t";
        }
        std::cout << "\n";
    }
    theta_file.close();

    // read matrix W
    matrix W;
    file_name = "W.csv";
    full_path = base_path + file_name;
    std::ifstream W_file(full_path);
    // Check that it has opened correctly
    if (!W_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open W.csv\n";
    }
    W = read_matrix(full_path, W);
    std::cout << "We print matrix W \n";
    for (auto i = 0; i <= 5; i ++) {
        for (auto j = 0; j <= 5; j++) {
            std::cout << W[i][j] << "\t";
        }
        std::cout << "\n";
    }
    W_file.close();

    // read matrix B_pred
    matrix B_pred;
    file_name = "Bpred.csv";
    full_path = base_path + file_name;
    std::ifstream B_pred_file(full_path);
    // Check that it has opened correctly
    if (!B_pred_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open Bpred.csv\n";
    }
    B_pred = read_matrix(full_path, B_pred);
    std::cout << "We print matrix B_pred \n";
    for (auto i = 0; i <= 5; i ++) {
        for (auto j = 0; j <= 5; j++) {
            std::cout << B_pred[i][j] << "\t";
        }
        std::cout << "\n";
    }
    B_pred_file.close();

    // read matrix thr1
    matrix thr1;
    file_name = "thr1.csv";
    full_path = base_path + file_name;
    std::ifstream thr1_file(full_path);
    // Check that it has opened correctly
    if (!thr1_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open thr1.csv\n";
    }
    thr1 = read_matrix(full_path, thr1);
    std::cout << "We print matrix thr1 \n";
    for (auto i = 0; i <= 5; i ++) {
        for (auto j = 0; j <= 5; j++) {
            std::cout << thr1[i][j] << "\t";
        }
        std::cout << "\n";
    }
    thr1_file.close();

    // read vector t_ij
    vector tij;
    file_name = "tij.csv";
    full_path = base_path + file_name;
    std::ifstream tij_file(full_path);
    // Check that it has opened correctly
    if (!tij_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open tij.csv\n";
    }
    tij = read_vector(full_path, tij);
    std::cout << "Printing tij" << "\n";
    for (auto i = 0; i <= 5; i ++) {
        std::cout << tij[i] << "\n";
    }
    tij_file.close();

    // read vector tg
    vector tg;
    file_name = "tg.csv";
    full_path = base_path + file_name;
    std::ifstream tg_file(full_path);

    // Check that it has opened correctly
    if (!tg_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open tg.csv\n";
    }
    tg = read_vector(full_path, tg);
    std::cout << "Printing tg" << "\n";
    for (auto i = 0; i <= 5; i ++) {
        std::cout << tg[i] << "\n";
    }
    tg_file.close();

    // read vector numeff
    vector numeff;
    file_name = "numeff.csv";
    full_path = base_path + file_name;
    std::ifstream numeff_file(full_path);

    // Check that it has opened correctly
    if (!numeff_file.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open numeff.csv\n";
    }
    numeff = read_vector(full_path, numeff);
    std::cout << "Printing numeff" << "\n";
    for (auto i = 0; i <= 5; i ++) {
        std::cout << numeff[i] << "\n";
    }
    numeff_file.close();

    // read int p
    int p;
    file_name = "p.csv";
    full_path = base_path + file_name;
    std::ifstream p_file(full_path);

    // Check that it has opened correctly
    if (!p_file.is_open()) {
        std::cerr << "Error opening file p.csv!" << std::endl;
    }

    p_file >> p;

    std::cout << "Printing p" << "\n";
    std::cout << p << "\n";

    // read int T
    int T;
    file_name = "T.csv";
    full_path = base_path + file_name;
    std::ifstream T_file(full_path);

    // Check that it has opened correctly
    if (!T_file.is_open()) {
        std::cerr << "Error opening file T.csv!" << std::endl;
    }

    T_file >> T;

    std::cout << "Printing T" << "\n";
    std::cout << T << "\n";

    T_file.close();

}



    /*
    std::ifstream tijfile("/Users/giuliadesanctis/circadian_genes/circadian_genes_bayesian/Data/tij.csv");
    //
    //std::vector<double> tg;

    if (!tijfile.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open tij.csv\n";
        return 1; // Exit with an error code
    }





    // Check that it has read tij correctly by printing first 20 values
    std::cout << "Printing t_ij" << "\n";
    for (auto i = 0; i <= 5; i ++) {
            std::cout << tij[i] << "\n";
    }


    std::ifstream tgfile("/Users/giuliadesanctis/circadian_genes/circadian_genes_bayesian/Data/tg.csv");

    if (!tgfile.is_open()) {
        // Check if the file opened successfully - this can be commented once the code is working well
        std::cerr << "Error: Could not open tg.csv\n";
        return 1; // Exit with an error code
    }

    while (std::getline(tgfile, line)) {
        try {
            tg.push_back(std::stod(line)); // Convert line to double and add to vector
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid line: '" << line << "' (not a valid number)\n";
        }
    }
    tgfile.close();

    std::cout << "Printing tg" << "\n";
    for (auto i = 0; i <= 5; i ++) {
        std::cout << tij[i] << "\n";
    }
}
     */





