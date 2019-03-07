//
//  main.cpp
//  405HW4
//
//  Created by Justin Tan on 2/5/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

int main(int argc, const char * argv[]) {
    //Q1
    cout << "Problem 1" << endl;
    double r1 = 0.05;
    double sigma1 = 0.24;
    double s1 = 32.0;
    double k1 = 30.0;
    double t = 6.0/12.0;
    
    vector<double> n{10,20,40,80,100,200,500};
    vector<double> c_a;
    vector<double> c_b;
    vector<double> c_c;
    vector<double> c_d;
    for(int i = 0; i < n.size(); i++){
        c_a.push_back(binomEuroCallA(r1, sigma1, s1, k1, t, n[i]));
        c_b.push_back(binomEuroCallB(r1, sigma1, s1, k1, t, n[i]));
        c_c.push_back(binomEuroCallJR(r1, sigma1, s1, k1, t, n[i]));
        c_d.push_back(binomEuroCallCRR(r1, sigma1, s1, k1, t, n[i]));
    }
    vector<vector<double>> data;
    data.push_back(n);
    data.push_back(c_a);
    data.push_back(c_b);
    data.push_back(c_c);
    data.push_back(c_d);
    toCSV(data, "q1.csv");
    cout << "All plots are in PDF." << endl;
    
    //Q2
    //price as of 11 am, Feb. 6
    cout << "\nProblem 2" << endl;
    double goog = 1115.23;
    double r2 = 0.02;
    double k2 = 1220;
    cout << "Stock data are from Yahoo Finance." << endl;
    vector<double> hist_p;
    ifstream file("GOOG.csv");
    string value = "";
    if(file.is_open()){
        getline(file, value, '\n');
        while(!file.eof()){
            for(int i = 0; i < 5; i++){
                getline(file, value, ',');
                value = "";
            }
            getline(file, value, ',');
            double price = atof(value.c_str());
            hist_p.push_back(price);
            getline(file, value, '\n');
            value = "";
        }
    }
    hist_p.pop_back();
    vector<double> hist_ret;
    for(int i = 1; i < hist_p.size(); i++){
        hist_ret.push_back(hist_p[i] / hist_p[i-1] - 1.0);
    }
    double sigma2 = calcStdDev(hist_ret)*sqrt(252);
    cout << "Estimated volatility: " << sigma2 << endl;
    double c2 = binomEuroCallCRR(r2, sigma2, goog, k2, 11.0/12.0, 500);
    cout << "Computed call option price: " << c2 << endl;
    //its within the bid-ask spread
    cout << "Actual call option price: 68.5" << endl;
    double sigma2_2 = 0.2375;
    cout << "Volatility to make call prices equal: " << sigma2_2 << endl;
    
    //Q3
    cout << "\nProblem 3\nAll plots are in PDF." << endl;
    double s3  = 49.0;
    double k3 = 50.0;
    double r3 = 0.03;
    double sigma3 = 0.2;
    double T3 = 0.3846;
    double mu3 = 0.14;
    vector<double> s_range;
    for(double i = 20.0; i<= 80; i += 2){
        s_range.push_back(i);
    }
    vector<double> t_range;
    for(double i = 0.01; i <= 0.3846; i +=0.01){
        t_range.push_back(i);
    }
    
    vector<double> delta_s;
    vector<double> delta_t;
    vector<double> thetaVal;
    vector<double> gammaVal;
    vector<double> vegaVal;
    vector<double> rhoVal;
    for(double i : s_range){
        delta_s.push_back(delta(r3, T3, sigma3, i, k3));
        thetaVal.push_back(theta(r3, T3, sigma3, i, k3));
        gammaVal.push_back(gamma(r3, T3, sigma3, i, k3));
        vegaVal.push_back(vega(r3, T3, sigma3, i, k3));
        rhoVal.push_back(theta(r3, T3, sigma3, i, k3));
    }
    for(double i : t_range){
        delta_t.push_back(delta(r3, i, sigma3, s3, k3));
    }
    vector<vector<double>> deltaData;
    deltaData.push_back(s_range);
    deltaData.push_back(t_range);
    deltaData.push_back(delta_s);
    deltaData.push_back(delta_t);
    deltaData.push_back(thetaVal);
    deltaData.push_back(gammaVal);
    deltaData.push_back(vegaVal);
    deltaData.push_back(rhoVal);
    toCSV(deltaData, "q3.csv");
    
    //Q4
    cout << "\nProblem 4\nAll plots are in PDF." << endl;
    double t4 = 1.0;
    double r4 = 0.05;
    double sigma4 = 0.3;
    double k4 = 100;
    vector<double> sRange4;
    for(double i = 80.0; i <= 120.0; i += 4.0){
        sRange4.push_back(i);
    }
    vector<double> euroPut;
    vector<double> amPut;
    for(double i : sRange4){
        euroPut.push_back(binomEuroPutCRR(r4, sigma4, i, k4, t4, 500));
        amPut.push_back(binomAmPutCRR(r4, sigma4, i, k4, t4, 500));
    }
    vector<vector<double>> dataQ4;
    dataQ4.push_back(sRange4);
    dataQ4.push_back(euroPut);
    dataQ4.push_back(amPut);
    toCSV(dataQ4, "q4.csv");
    
    //Q5
    cout << "\nProblem 5\nAll plots are in PDF." << endl;
    double r5 = 0.05;
    double sigma5 = 0.24;
    double s5 = 32.0;
    double k5 = 30.0;
    double t5 = 0.5;
    vector<int> n5{10,15,20,40,70,80,100,200,500};
    vector<double> n5Double{10,15,20,40,70,80,100,200,500};
    vector<double> trinom;
    vector<double> trinomLog;
    for(int i : n5){
        trinom.push_back(trinomEuroCall(r5, sigma5, s5, k5, t5, i));
        trinomLog.push_back(trinomEuroCallLog(r5, sigma5, s5, k5, t5, i));
    }
    vector<vector<double>> dataQ5;
    dataQ5.push_back(n5Double);
    dataQ5.push_back(trinom);
    dataQ5.push_back(trinomLog);
    toCSV(dataQ5, "q5.csv");
    
    //Q6
    cout << "\nProblem 6." << endl;
    double s6 = 50;
    double k6 = 55;
    double T6 = 2;
    double r6 = 0.02;
    double sigma6 = 0.2;
    int n6 = 1000;
    int base1 = 2;
    int base2 = 5;
    cout << "Price : " << euroCall(r6, T6, sigma6, s6, k6, n6, base1, base2) << endl;
    
    return 0;
    

}
