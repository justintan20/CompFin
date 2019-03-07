//
//  main.cpp
//  405HW9
//
//  Created by Justin Tan on 3/6/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

int main(int argc, const char * argv[]) {
    cout << "Problem 1" << endl;
    double loan = 100000;
    double wac = 0.08;
    unsigned int T_years = 30;
    double r0 = 0.078;
    double kappa = 0.6;
    double r_bar = 0.08;
    double sigma = 0.12;
    int m = 10000;
    //a
    double price = priceMBSNumerix(T_years, loan, wac, r0, kappa, r_bar, sigma, m);
    
    cout << price << endl;
    
    
    return 0;
}
