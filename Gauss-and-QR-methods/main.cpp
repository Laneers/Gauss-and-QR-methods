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

    bool isSingular = false;

    Gauss_method<float>(filename, &isSingular);
    Gauss_method<double>(filename, &isSingular);

    QR_method<float>(filename, &isSingular);
    QR_method<double>(filename, &isSingular);

    ConditionNumber<float>(filename, "Gauss", &isSingular);
    ConditionNumber<double>(filename, "Gauss", &isSingular);

    ConditionNumber<float>(filename, "QR", &isSingular);
    ConditionNumber<double>(filename, "QR", &isSingular);

    std::cout << "\nAll methods completed successfully.\n";
    std::cout << "\nTesting\n";


    return 0;
}