#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include "Norms.h"

template <typename T>
void QR_method(const std::string filename, bool* isSingular, bool perturbed = false) {
    //Чтение данных из файла
    std::ifstream file(filename + ".txt");
    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return;
    }

    int n;
    file >> n;
    T** A_original = new T * [n + 1];
    T** R = new T * [n + 1];
    T** Q = new T * [n + 1];
    T** T_mat = new T * [n + 1];
    for (int i = 1; i <= n; i++) {
        A_original[i] = new T[n + 1];
        R[i] = new T[n + 1];
        Q[i] = new T[n + 1];
        T_mat[i] = new T[n + 1];
    }
    T* b = new T[n + 1];

    T val;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n + 1; j++) {
            file >> val;
            if (j <= n) {
                A_original[i][j] = val;
                R[i][j] = val;
            }
            if (j == n + 1) {
                b[i] = val;
            }
        }
    }
    file.close();

    //Инициализация Q как единичной матрицы
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 1; i < n; i++) {
        for (int j = 1; j <= n; j++) {
            T_mat[i][j] = (i == j) ? 1.0 : 0.0;
        }
        for (int j = i + 1; j <= n; j++) {
            if (std::abs(R[j][i]) < std::numeric_limits<T>::epsilon() * 10) {
                continue;
            }
            T denominator = std::sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
            T c = R[i][i] / denominator;
            T s = R[j][i] / denominator;

            for (int k = i; k <= n; k++) {
                T r_ik = R[i][k];
                T r_jk = R[j][k];
                R[i][k] = c * r_ik + s * r_jk;
                R[j][k] = -s * r_ik + c * r_jk;
            }

            for (int k = 1; k <= n; k++) {
                T q_ki = Q[k][i];
                T q_kj = Q[k][j];
                Q[k][i] = c * q_ki + s * q_kj;
                Q[k][j] = -s * q_ki + c * q_kj;
            }
        }
    }

    T norm_A = matrix_norm_inf(A_original, n);
    T threshold = std::numeric_limits<T>::epsilon() * norm_A * n;
    /*T threshold = typeid(T).name() == "double" ? 1e-12 : 1e-6;*/
    for (int i = 1; i <= n; i++) {
        if (std::fabs(R[i][i]) < threshold) {
            *isSingular = true;
            break;
        }
    }

    std::string filename_ans = std::string("QR_answer") + filename[4] + '_' + typeid(T).name() + ".txt";
    if (perturbed) {
        filename_ans = std::string("QR_answer") + filename[4] + "_perturbed_" + typeid(T).name() + ".txt";
    }

    if (*isSingular) {
        std::cout << "\nErrorQ: matrix in " + filename + " is degenerate\n";
        std::ofstream out(filename_ans);
        out << "Error: matrix is degenerate";
        out.close();

        for (int i = 1; i <= n; i++) {
            delete[] A_original[i];
            delete[] R[i];
            delete[] Q[i];
            delete[] T_mat[i];
        }
        delete[] A_original;
        delete[] R;
        delete[] Q;
        delete[] T_mat;
        delete[] b;
        return;
    }

    T* y = new T[n + 1];
    for (int i = 1; i <= n; i++) {
        y[i] = 0.0;
        for (int j = 1; j <= n; j++) {
            y[i] += Q[j][i] * b[j];     //транспонированная Q
        }
    }
    T* x = new T[n + 1];
    x[n] = y[n] / R[n][n];
    for (int i = n - 1; i >= 1; i--) {
        T sum = 0.0;
        for (int j = i + 1; j <= n; j++) {
            sum += R[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / R[i][i];
    }

    //Невязка
    T* b1 = new T[n + 1];
    for (int i = 1; i <= n; i++) {
        b1[i] = 0.0;
        for (int j = 1; j <= n; j++) {
            b1[i] += A_original[i][j] * x[j];
        }
    }
    T* err = new T[n + 1];
    for (int i = 1; i <= n; i++) {
        err[i] = b[i] - b1[i];
    }
    double normL2 = 0.0;
    for (int i = 1; i <= n; i++) {
        normL2 += err[i] * err[i];
    }
    normL2 = sqrt(normL2);

    //Выводим решения
    std::ofstream out(filename_ans);
    for (int i = 1; i <= n; i++) {
        out << x[i] << ' ';
    }
    out << "\nNormL2 = " << normL2;
    out << "\n";
    out << "\nMatrix Q\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            out << Q[i][j] << " ";
        }
        out << "\n";
    }
    out << "\nMatrix R\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            out << R[i][j] << " ";
        }
        out << "\n";
    }
    out.close();

    for (int i = 1; i <= n; i++) {
        delete[] A_original[i];
        delete[] R[i];
        delete[] Q[i];
        delete[] T_mat[i];
    }
    delete[] A_original;
    delete[] R;
    delete[] Q;
    delete[] T_mat;
    delete[] b;
    delete[] y;
    delete[] x;
    delete[] b1;
    delete[] err;
}