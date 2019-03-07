//
//  main.cpp
//  405HW6
//
//  Created by Justin Tan on 2/20/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

int main(int argc, const char * argv[]) {
    cout << "Problem 1" << endl;
    cout << "Running...";
    double t = 1.0;
    double r = 0.03;
    double s0 = 98;
    double k = 100;
    int m = 10000;
    int n = 100;
    vector<double> sigmas;
    for(double i = 0.12; i <= 0.48; i += 0.04){
        sigmas.push_back(i);
    }
    vector<double> calls;
    vector<double> puts;
    for(double vol : sigmas){
        calls.push_back(fixedStrikeLookbackCall(r, s0, k, vol, t, m, n));
        puts.push_back(fixedStrikeLookbackPut(r, s0, k, vol, t, m, n));
        cout << ".";
    }
    vector<vector<double>> q1Data{sigmas, calls, puts};
    toCSV(q1Data, "q1.csv");
    cout << "\n Data is outputted as \"q1.csv\" and all plots are in PDF." << endl;
    
    cout << "\nProblem 2" << endl;
    cout << "Running...";
    double value = 0;
    double prob = 0;
    double Et = 0;
    vector<vector<double>> values_lambda1_const;
    vector<vector<double>> probs_lambda1_const;
    vector<vector<double>> Ets_lambda1_const;
    for(double lambda2 = 0; lambda2 <= 0.8; lambda2 += 0.1){
        vector<double> values;
        vector<double> probs;
        vector<double> Ets;
        for(int T = 3; T <= 8; T++){
            Proj6_2function(0.2, lambda2, T, value, prob, Et);
            cout << ".";
            values.push_back(value);
            probs.push_back(prob);
            Ets.push_back(Et);
        }
        values_lambda1_const.push_back(values);
        probs_lambda1_const.push_back(probs);
        Ets_lambda1_const.push_back(Ets);
    }


    vector<vector<double>> values_lambda2_const;
    vector<vector<double>> probs_lambda2_const;
    vector<vector<double>> Ets_lambda2_const;
    for(double lambda1 = 0.05; lambda1 <= 0.4; lambda1 += 0.05){
        vector<double> values;
        vector<double> probs;
        vector<double> Ets;
        for(int T = 3; T <= 8; T++){
            cout << ".";
            Proj6_2function(lambda1, 0.4, T, value, prob, Et);
            values.push_back(value);
            probs.push_back(prob);
            Ets.push_back(Et);
        }
        values_lambda2_const.push_back(values);
        probs_lambda2_const.push_back(probs);
        Ets_lambda2_const.push_back(Ets);
    }

    toCSV(values_lambda1_const, "lambda1_val.csv");
    toCSV(probs_lambda1_const, "lambda1_prob.csv");
    toCSV(Ets_lambda1_const, "lambda1_et.csv");
    toCSV(values_lambda2_const, "lambda2_val.csv");
    toCSV(probs_lambda2_const, "lambda2_prob.csv");
    toCSV(Ets_lambda2_const, "lambda2_et.csv");
    
    cout << "\n All data outputted as CSV files, all plots in PDF." << endl;
    cout << "Default Parameter values:" << endl;
    Proj6_2function(0.2, 0.4, 5, value, prob, Et);
    cout << "Option Value: " << value << endl;
    cout << "Default Probability: " << prob << endl;
    cout << "Expected Exercise Time: " << Et << endl;

    return 0;
}
