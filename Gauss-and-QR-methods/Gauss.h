#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include "Norms.h"

template <typename T>
void Gauss_method(const std::string filename, bool* isSingular, bool perturbed = false) {
    //Чтение данных из файла
    std::ifstream file(filename + ".txt");

    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return;
    }

    int n;
    file >> n;
    int m = n + 1;
    T** A_original = new T * [n + 1];
    T** A = new T * [n + 1];
    for (int i = 1; i <= n; i++) {
        A_original[i] = new T[n + 1];
        A[i] = new T[m + 1];    //расширенная матрица
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

    //Выбор главного элемента
    int index;
    for (int i = 1; i < n; i++) {
        T maxEl = fabs(A[i][i]);
        index = i;
        for (int j = i + 1; j <= n; j++) {
            if (fabs(A[j][i]) > maxEl) {
                maxEl = fabs(A[j][i]);
                index = j;
            }
        }
        for (int k = 1; k <= m; k++) {     //меняем строки
            T temp = A[i][k];
            A[i][k] = A[index][k];
            A[index][k] = temp;
        }

        //Прямой ход
        T tmp;
        for (int l = i + 1; l <= n; l++) {
            tmp = A[l][i] / A[i][i];
            for (int t = i; t <= m; t++) {
                A[l][t] -= tmp * A[i][t];
            }
        }
    }

    T norm_A = matrix_norm_inf(A_original, n);
    T threshold = std::numeric_limits<T>::epsilon() * norm_A * n;
    /*T threshold = typeid(T).name() == "double" ? 1e-12 : 1e-6;*/
    for (int i = 1; i <= n; i++) {
        if (std::fabs(A[i][i]) < threshold) {
            *isSingular = true;
            break;
        }
    }

    std::string filename_ans = std::string("Gauss_answer") + filename[4] + '_' + typeid(T).name() + ".txt";
    if(perturbed) {
        filename_ans = std::string("Gauss_answer") + filename[4] + "_perturbed_" + typeid(T).name() + ".txt";
	}

    if (*isSingular) {
        std::cout << "\nErrorG: matrix in " + filename + " is degenerate\n";
        std::ofstream out(filename_ans);
        out << "Error: matrix is degenerate";
        out.close();

        for (int i = 1; i <= n; i++) {
            delete[] A[i];
            delete[] A_original[i];
        }
        delete[] A;
        delete[] A_original;
        delete[] b;
        return;
    }

    //Обратный ход
    T* x = new T[n + 1];
    x[n] = A[n][m] / A[n][n];
    for (int i = n - 1; i >= 1; i--) {
        x[i] = A[i][m];
        for (int j = i + 1; j <= n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] = x[i] / A[i][i];
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
    out.close();

    for (int i = 1; i <= n; i++) {
        delete[] A[i];
        delete[] A_original[i];
    }
    delete[] A;
    delete[] A_original;
    delete[] x;
    delete[] b;
    delete[] b1;
    delete[] err;
}