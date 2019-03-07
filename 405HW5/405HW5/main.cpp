//
//  main.cpp
//  405HW5
//
//  Created by Justin Tan on 2/14/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

int main(int argc, const char * argv[]) {
    cout << "Problem 1" << endl;
    cout << "Running...";
    double sigma = 0.2;
    double r = 0.06;
    int m = 100000;
    int n = 100;
    double k = 40.0;
    vector<double> time{0.5, 1.0, 2.0};
    vector<vector<double>> results;
    for(int i = 0; i < 3; i++){
        vector<double> method_result;
        switch (i) {
            case 0:
                for(int num = 2; num <= 4; num++){
                    for(double price = 36.0; price <= 44.0; price += 4.0){
                        for(double t : time){
                            method_result.push_back(LSMC(price, k, r, sigma, t, m, n, num, "hermite"));
                            cout << ".";
                        }
                    }
                }
                break;
            case 1:
                for(int num = 2; num <= 4; num++){
                    for(double price = 36.0; price <= 44.0; price += 4.0){
                        for(double t : time){
                            method_result.push_back(LSMC(price, k, r, sigma, t, m, n, num, "laguerre"));
                            cout << ".";
                        }
                    }
                }
                break;
            case 2:
                for(int num = 2; num <= 4; num++){
                    for(double price = 36.0; price <= 44.0; price += 4.0){
                        for(double t : time){
                            method_result.push_back(LSMC(price, k, r, sigma, t, m, n, num, "monomial"));
                            cout << ".";
                        }
                    }
                }
                break;
        }
        results.push_back(method_result);
    }
    toCSV(results, "put_prices.csv");
    cout << "\nAll results are in file \"put_prices.csv\"." << endl;
    
    cout << "\nProblem 2" << endl;
    double s0_2 = 65.0;
    double sigma_2 = 0.2;
    double r_2 = 0.06;
    double t_2 = 0.2;
    double T_2 = 1;
    double price_euro_forward = forwardStartEuroPut(s0_2, sigma_2, r_2, t_2, T_2, 100, 10000);
    cout << "Price of Forward-Start European Put: " << price_euro_forward << endl;
    double price_american_forward = forwardStartAmericanPut(s0_2, sigma_2, r_2, T_2, t_2, 10000, 100, 3, "monomial");
    cout << "Price of Forward-Start American Put (Using Least-Square Monte Carlo): " << price_american_forward << endl;
    return 0;
}
