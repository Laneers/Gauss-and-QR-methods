#include <iostream>
#include <fstream>
#include <cmath>

#include "Gauss.h"
#include "QR.h"
#include "ConditionNumber.h"

//double EPS_fl = 1e-6;
//double EPS_db = 3*1e-16;

int main() {

    std::string filename;
    std::cin >> filename;

    bool isSingular_double = false;
    bool isSingular_float = false;

    Gauss_method<float>(filename, &isSingular_float);
    QR_method<float>(filename, &isSingular_float);
    ConditionNumber<float>(filename, "Gauss", &isSingular_float);
    ConditionNumber<float>(filename, "QR", &isSingular_float);

    Gauss_method<double>(filename, &isSingular_double);
    QR_method<double>(filename, &isSingular_double);
    ConditionNumber<double>(filename, "Gauss", &isSingular_double);
    ConditionNumber<double>(filename, "QR", &isSingular_double);

    std::cout << "\nAll methods completed successfully.\n";

    return 0;
}