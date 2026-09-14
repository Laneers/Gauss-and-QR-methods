#pragma once

#include <cmath>

template <typename T>
T vector_norm_l1(T* v, int n) {
    T sum = 0;
    for (int i = 1; i <= n; i++) {
        sum += std::fabs(v[i]);
    }
    return sum;
}

template <typename T>
T matrix_norm_l1(T** A, int n) {
    T max_col_sum = 0;
    for (int j = 1; j <= n; j++) {      //идем по столбцам
        T col_sum = 0;
        for (int i = 1; i <= n; i++) {      //cуммируем строки в столбце
            col_sum += std::fabs(A[i][j]);
        }
        if (col_sum > max_col_sum) {
            max_col_sum = col_sum;
        }
    }
    return max_col_sum;
}

template <typename T>
T vector_norm_inf(T* v, int n) {
    T max_val = 0;
    for (int i = 1; i <= n; i++) {
        if (std::fabs(v[i]) > max_val) {
            max_val = std::fabs(v[i]);
        }
    }
    return max_val;
}

template <typename T>
T matrix_norm_inf(T** A, int n) {
    T max_row_sum = 0;
    for (int i = 1; i <= n; i++) {
        T row_sum = 0;
        for (int j = 1; j <= n; j++) {
            row_sum += std::fabs(A[i][j]);
        }
        if (row_sum > max_row_sum) {
            max_row_sum = row_sum;
        }
    }
    return max_row_sum;
}