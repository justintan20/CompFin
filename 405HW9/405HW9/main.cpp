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

    cout << "Price of Mortgage: " << price << endl;

    vector<double> kappa_range;
    vector<double> r_bar_range;
    vector<double> sigma_range;
    for(double k = 0.3; k <= 0.9; k += 0.1){
        kappa_range.push_back(k);
        r_bar_range.push_back(k/10.0);
    }
    for(double vol = 0.1; vol <= 0.2; vol += 0.01){
        sigma_range.push_back(vol);
    }
    vector<double> prices_kappa;
    for(double k : kappa_range){
        prices_kappa.push_back(priceMBSNumerix(T_years, loan, wac, r0, k, r_bar, sigma, m));
    }
    vector<double> prices_rbar;
    for(double r : r_bar_range){
        prices_rbar.push_back(priceMBSNumerix(T_years, loan, wac, r0, kappa, r, sigma, m));
    }
    vector<double> prices_sigma;
    for(double s : sigma_range){
        prices_sigma.push_back(priceMBSNumerix(T_years, loan, wac, r0, kappa, r_bar, s, m));
    }
    vector<vector<double>> kappa_data{kappa_range,prices_kappa};
    vector<vector<double>> rbar_data{r_bar_range,prices_rbar};
    vector<vector<double>> sigma_data{sigma_range,prices_sigma};
    toCSV(kappa_data, "kappa.csv");
    toCSV(rbar_data, "rbar.csv");
    toCSV(sigma_data, "sigma.csv");
    cout << "All data are outputed as CSV files, all plots are in PDF." << endl;
    
    
    cout << "\nProblem 2" << endl;
    double marketPrice = 110000;
    double x = oas(T_years, loan, wac, r0, kappa, r_bar, sigma, m, marketPrice);
    cout << "Option Adjusted Spread: " << x << endl;
    
    cout << "\nProblem 3" << endl;
    double y = 0.0005;
    double p_plus = priceMBSNumerixWithOAS(T_years, loan, wac, r0, kappa, r_bar, sigma, m, x+y);
    double p_minus = priceMBSNumerixWithOAS(T_years, loan, wac, r0, kappa, r_bar, sigma, m, x-y);
    double duration = (p_minus - p_plus)/(2.0*y*marketPrice);
    double convexity = (p_plus+p_minus-2.0*marketPrice)/(2.0*marketPrice*y*y);
    cout << "OAS-Duration: " << duration << endl;
    cout << "OAS-Convexity: " << convexity << endl;
    
    
    return 0;
}
