#include <iostream>
#include <fstream>
#include <cmath>

#include "Gauss.h"
#include "QR.h"
#include "ConditionNumber.h"

//double EPS_fl = 1e-6;
//double EPS_db = 3*1e-16;

void generate_perturbed_file(const std::string filename, double* delta_b) {
    //Чтение данных из файла
    std::ifstream file(filename + ".txt");
    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return;
    }

    int n;
    file >> n;
    double** A_original = new double * [n + 1];
    for (int i = 1; i <= n; i++) {
        A_original[i] = new double[n + 1];
    }
    double* b = new double[n + 1];

    double val;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n + 1; j++) {
            file >> val;
            if (j <= n) {
                A_original[i][j] = val;
            }
            if (j == n + 1) {
                b[i] = val;
            }
        }
    }
    file.close();

    double* b_perturbed = new double[n + 1];
    for (int i = 1; i <= n; i++) {
        double sign = (std::rand() % 2 == 0) ? 1.0 : -1.0;
        delta_b[i] = sign * 0.01;
        b_perturbed[i] = b[i] + delta_b[i];
    }
    std::string perturbed_filename = filename + "_perturbed";
    std::ofstream file_p(perturbed_filename + ".txt");
    file_p << n << "\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            file_p << A_original[i][j] << " ";
        }
        file_p << b_perturbed[i] << "\n";
    }
    file_p.close();

    for (int i = 1; i <= n; i++) {
        delete[] A_original[i];
    }
    delete[] A_original;
    delete[] b;
    delete[] b_perturbed;
}

int main() {

    std::string filename;
    std::cin >> filename;

    //Чтение данных из файла
    std::ifstream file(filename + ".txt");
    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return 1;
    }
    int n;
    file >> n;

    bool isSingular_double = false;
    bool isSingular_float = false;
    bool isPerturbedSingular = false;

    double* delta_b_double = new double[n + 1];
    generate_perturbed_file(filename, delta_b_double); 

    // Создаем float-копию возмущения
    float* delta_b_float = new float[n + 1];
    for (int i = 1; i <= n; i++) {
        delta_b_float[i] = static_cast<float>(delta_b_double[i]);
    }

    Gauss_method<float>(filename, &isSingular_float);
    QR_method<float>(filename, &isSingular_float);
    ConditionNumber<float>(filename, "Gauss", delta_b_float, &isPerturbedSingular);
    isPerturbedSingular = false;
    ConditionNumber<float>(filename, "QR", delta_b_float, &isPerturbedSingular);

    Gauss_method<double>(filename, &isSingular_double);
    QR_method<double>(filename, &isSingular_double);
    isPerturbedSingular = false;
    ConditionNumber<double>(filename, "Gauss", delta_b_double, &isPerturbedSingular);
    isPerturbedSingular = false;
    ConditionNumber<double>(filename, "QR", delta_b_double, &isPerturbedSingular);

    return 0;
}