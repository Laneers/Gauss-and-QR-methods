#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include "Norms.h"

template <typename T>
T** invert_matrix(const std::string& filename, bool* isSingular) {
    //Чтение данных из файла
    std::ifstream file(filename + ".txt");
    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return nullptr;
    }

    int n;
    file >> n;
    int m = n + 1;
    T** A_original = new T * [n + 1];
    T** A_inverse = new T * [n + 1];
    T** A = new T * [n + 1];
    for (int i = 1; i <= n; i++) {
        A_original[i] = new T[n + 1];
        A_inverse[i] = new T[n + 1];
        A[i] = new T[m + 1];
    }
    T* b = new T[n + 1];

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            file >> A[i][j];
            if (j <= n) {
                A_original[i][j] = A[i][j];
            }
            if (j == m) {
                b[i] = A[i][j];
            }
        }
    }
    file.close();

    T norm_A = matrix_norm_inf(A_original, n); 
    T threshold = std::numeric_limits<T>::epsilon() * norm_A * n;

    for (int col = 1; col <= n; col++) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                A[i][j] = A_original[i][j];
            }
            A[i][m] = (i == col) ? 1.0 : 0.0;
        }

        int index;
        for (int i = 1; i < n; i++) {
            T maxEl = std::fabs(A[i][i]);
            index = i;
            for (int j = i + 1; j <= n; j++) {
                if (std::fabs(A[j][i]) > maxEl) {
                    maxEl = std::fabs(A[j][i]);
                    index = j;
                }
            }
            for (int k = 1; k <= m; k++) {     // Меняем строки
                T temp = A[i][k];
                A[i][k] = A[index][k];
                A[index][k] = temp;
            }

            T tmp;
            for (int l = i + 1; l <= n; l++) {
                tmp = A[l][i] / A[i][i];
                for (int t = i; t <= m; t++) {
                    A[l][t] -= tmp * A[i][t];
                }
            }
        }

        if (col == 1) {
            for (int i = 1; i <= n; i++) {
                if (std::fabs(A[i][i]) < threshold) {
                    *isSingular = true;
                    break;
                }
            }
            if (*isSingular) {
                std::cout << "\nErrorI: matrix in " + filename + " is degenerate\n";
                for (int i = 1; i <= n; i++) {
                    delete[] A[i];
                    delete[] A_original[i];
                    delete[] A_inverse[i];
                }
                delete[] A;
                delete[] A_inverse;
                delete[] A_original;
                delete[] b;
                return nullptr;
            }
        }

        T* x = new T[n + 1];
        x[n] = A[n][m] / A[n][n];
        for (int i = n - 1; i >= 1; i--) {
            x[i] = A[i][m];
            for (int j = i + 1; j <= n; j++) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] = x[i] / A[i][i];
        }

        for (int i = 1; i <= n; i++) {
            A_inverse[i][col] = x[i];
        }

        delete[] x;
    }

    for (int i = 1; i <= n; i++) {
        delete[] A[i];
    }
    delete[] A;

    return A_inverse;
}

template <typename T>
void ConditionNumber(const std::string filename, const std::string method_name, T* delta_b, bool* isPerturbedSingular) {
    //Чтение данных из файла
    std::ifstream file(filename + ".txt");
    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return;
    }

    int n;
    file >> n;
    T** A_original = new T * [n + 1];
    for (int i = 1; i <= n; i++) {
        A_original[i] = new T[n + 1];
    }
    T* b = new T[n + 1];

    T val;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n + 1; j++) {
            file >> val;
            if (j <= n) {
                A_original[i][j] = val;
            }
            else {
                b[i] = val;
            }
        }
    }
    file.close();

    //Подсчёт норм
    bool isInvertSingular = false;
    T** A_inv = invert_matrix<T>(filename, &isInvertSingular);

    T cond_true_l1 = -1;
    T cond_true_inf = -1;

    if (!(isInvertSingular) && A_inv != nullptr) {
        T norm_A_l1 = matrix_norm_l1(A_original, n);
        T norm_A_inv_l1 = matrix_norm_l1(A_inv, n);
        cond_true_l1 = norm_A_l1 * norm_A_inv_l1;

        T norm_A_inf = matrix_norm_inf(A_original, n);
        T norm_A_inv_inf = matrix_norm_inf(A_inv, n);
        cond_true_inf = norm_A_inf * norm_A_inv_inf;

        for (int i = 1; i <= n; i++) {
            delete[] A_inv[i];
        }
        delete[] A_inv;
    }

    if (method_name == "Gauss") {
        Gauss_method<T>(filename + "_perturbed", isPerturbedSingular, true);
    }
    if (method_name == "QR") {
        QR_method<T>(filename + "_perturbed", isPerturbedSingular, true);
    }

    T* x_perturbed = new T[n + 1];
    std::string filename_pertubed_ans = std::string(method_name) + std::string("_answer") + filename[4] + "_perturbed_" + typeid(T).name() + ".txt";
    std::ifstream ans_file_pertubed(filename_pertubed_ans);

    if (!ans_file_pertubed.is_open()) {
        std::cerr << "Error: failed to open perturbed answer file\n";
        delete[] x_perturbed;
        for (int i = 1; i <= n; i++) {
            delete[] A_original[i];
        }
        delete[] A_original;
        delete[] b;
        return;
    }

    //Проверяем, не выродилась ли матрица при расчете
    if (*isPerturbedSingular) {
        std::cout << "\nErrorC: pertubed matrix is degenerate\n";
        delete[] x_perturbed;
        for (int i = 1; i <= n; i++) {
            delete[] A_original[i];
        }
        delete[] A_original;
        delete[] b;
        return;
    }
    
    for (int i = 1; i <= n; i++) {
        ans_file_pertubed >> x_perturbed[i];
    }
    ans_file_pertubed.close();

    std::ifstream ans_file(std::string(method_name) + "_answer" + filename[4] + '_' + typeid(T).name() + ".txt");

    if (!ans_file.is_open()) {
        std::cerr << "Error: failed to open answer file\n";
        delete[] x_perturbed;
        for (int i = 1; i <= n; i++) {
            delete[] A_original[i];
        }
        delete[] A_original;
        delete[] b;
        return;
    }

    T* x_original = new T[n + 1];
    for (int i = 1; i <= n; i++) {
        ans_file >> x_original[i];
    }
    ans_file.close();

    //Вычисляем абсолютные погрешности
    T* delta_x = new T[n + 1];
    for (int i = 1; i <= n; i++) {
        delta_x[i] = x_perturbed[i] - x_original[i];
    }

    //Вычисляем нормы для анализа
    T norm_b_l1 = vector_norm_l1(b, n);
    T norm_delta_b_l1 = vector_norm_l1(delta_b, n);
    T norm_x_l1 = vector_norm_l1(x_original, n);
    T norm_delta_x_l1 = vector_norm_l1(delta_x, n);
    T norm_b_inf = vector_norm_inf(b, n);
    T norm_delta_b_inf = vector_norm_inf(delta_b, n);
    T norm_x_inf = vector_norm_inf(x_original, n);
    T norm_delta_x_inf = vector_norm_inf(delta_x, n);

    T rel_err_b_l1 = norm_delta_b_l1 / norm_b_l1;
    T rel_err_x_l1 = norm_delta_x_l1 / norm_x_l1;
    T rel_err_b_inf = norm_delta_b_inf / norm_b_inf;
    T rel_err_x_inf = norm_delta_x_inf / norm_x_inf;

    //Оценка коэффициента усиления ошибки
    T cond_est_l1 = rel_err_x_l1 / rel_err_b_l1;
    T cond_est_inf = rel_err_x_inf / rel_err_b_inf;

    std::cout << "\nCondition analysis of " << method_name << " for " << typeid(T).name() << " in L1 " << "\n";
    if (cond_true_l1 != -1) {
        std::cout << "True cond(A)_l1: " << cond_true_l1 << "\n";
    }
    else {
        std::cout << "Condition cannot be compute\n";
    }
    std::cout << "Estimation cond(A)_l1: " << cond_est_l1 << "\n";
    std::cout << "Relative error b: " << (norm_delta_b_l1 / norm_b_l1) << "\n";
    std::cout << "Relative error x: " << (norm_delta_x_l1 / norm_x_l1) << "\n";

    std::cout << "\nCondition analysis of " << method_name << " for " << typeid(T).name() << " in inf " << "\n";
    if (cond_true_inf != -1) {
        std::cout << "True cond(A)_inf: " << cond_true_inf << "\n";
    }
    else {
        std::cout << "Condition cannot be compute\n";
    }
    std::cout << "Estimation cond(A)_inf: " << cond_est_inf << "\n";
    std::cout << "Relative error b: " << (norm_delta_b_inf / norm_b_inf) << "\n";
    std::cout << "Relative error x: " << (norm_delta_x_inf / norm_x_inf) << "\n";

    for (int i = 1; i <= n; i++) {
        delete[] A_original[i];
    }
    delete[] A_original;
    delete[] b;
    delete[] x_perturbed;
    delete[] delta_x;
}