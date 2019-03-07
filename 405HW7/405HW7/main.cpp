//
//  main.cpp
//  405HW7
//
//  Created by Justin Tan on 2/26/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

int main(int argc, const char * argv[]) {
    cout << "Problem 1: European Put Prices" << endl;
    double s0 = 10;
    double s_low = 4;
    double s_high = 16;
    double sigma = 0.2;
    double r = 0.04;
    double delta_t = 0.002;
    double T = 0.5;
    double k = 10;
    double delta_x1 = sigma * sqrt(delta_t);
    double delta_x2 = sigma * sqrt(3*delta_t);
    double delta_x3 = sigma * sqrt(4*delta_t);
    cout << "Explicit Finite Difference:" << endl;
    double Pa_1 = euroPutEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x1, k);
    cout << "Pa_1: " << Pa_1 << endl;
    double Pa_2 = euroPutEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x2, k);
    cout << "Pa_2: " << Pa_2 << endl;
    double Pa_3 = euroPutEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x3, k);
    cout << "Pa_3: " << Pa_3 << endl;
    cout << "Implicit Finite Difference:" << endl;
    double Pb_1 = euroPutIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x1, k);
    cout << "Pb_1: " << Pb_1 << endl;
    double Pb_2 = euroPutIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x2, k);
    cout << "Pb_2: " << Pb_2 << endl;
    double Pb_3 = euroPutIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x3, k);
    cout << "Pb_3: " << Pb_3 << endl;
    cout << "Crank-Nicolson Finite Difference:" << endl;
    double Pc_1 = euroPutCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x1, k);
    cout << "Pc_1: " << Pc_1 << endl;
    double Pc_2 = euroPutCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x2, k);
    cout << "Pc_2: " << Pc_2 << endl;
    double Pc_3 = euroPutCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_x3, k);
    cout << "Pc_3: " << Pc_3 << endl;
    cout << "Black-Scholes Formula:" << endl;
    double P_BS = euroPutBS(r, T, sigma, s0, k);
    cout << "P: " << P_BS << endl;
    
    cout << "\nProblem 2: American Option Prices" << endl;
    double delta_s1 = 0.25;
    double delta_s2 = 1;
    double delta_s3 = 1.25;
    cout << "Put Options:" << endl;
    cout << "Explicit Finite Difference:" << endl;
    double Pa1 = amPutEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s1, k);
    cout << "Pa_1: " << Pa1 << endl;
    double Pa2 = amPutEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s2, k);
    cout << "Pa_2: " << Pa2 << endl;
    double Pa3 = amPutEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s3, k);
    cout << "Pa_3: " << Pa3 << endl;
    cout << "Implicit Finite Difference:" << endl;
    double Pb1 = amPutIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s1, k);
    cout << "Pb_1: " << Pb1 << endl;
    double Pb2 = amPutIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s2, k);
    cout << "Pb_2: " << Pb2 << endl;
    double Pb3 = amPutIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s3, k);
    cout << "Pb_3: " << Pb3 << endl;
    cout << "Crank-Nicolson Finite Difference:" << endl;
    double Pc1 = amPutCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s1, k);
    cout << "Pc_1: " << Pc1 << endl;
    double Pc2 = amPutCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s2, k);
    cout << "Pc_2: " << Pc2 << endl;
    double Pc3 = amPutCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s3, k);
    cout << "Pc_3: " << Pc3 << endl;
    cout << "Call Options:" << endl;
    cout << "Explicit Finite Difference:" << endl;
    double Ca1 = amCallEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s1, k);
    cout << "Ca_1: " << Ca1 << endl;
    double Ca2 = amCallEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s2, k);
    cout << "Ca_2: " << Ca2 << endl;
    double Ca3 = amCallEFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s3, k);
    cout << "Ca_3: " << Ca3 << endl;
    cout << "Implicit Finite Difference:" << endl;
    double Cb1 = amCallIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s1, k);
    cout << "Cb_1: " << Cb1 << endl;
    double Cb2 = amCallIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s2, k);
    cout << "Cb_2: " << Cb2 << endl;
    double Cb3 = amCallIFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s3, k);
    cout << "Cb_3: " << Cb3 << endl;
    cout << "Crank-Nicolson Finite Difference:" << endl;
    double Cc1 = amCallCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s1, k);
    cout << "Cc_1: " << Cc1 << endl;
    double Cc2 = amCallCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s2, k);
    cout << "Cc_2: " << Cc2 << endl;
    double Cc3 = amCallCNFD(s0, s_low,s_high,r, sigma, T, delta_t, delta_s3, k);
    cout << "Cc_3: " << Cc3 << endl;
    
    vector<double> call_efd;
    vector<double> call_ifd;
    vector<double> call_cnfd;
    for(double i = s_low; i <= s_high; i+= delta_s1){
        call_efd.push_back(amCallEFD(i, s_low,s_high,r, sigma, T, delta_t, delta_s1, k));
        call_ifd.push_back(amCallIFD(i, s_low,s_high,r, sigma, T, delta_t, delta_s1, k));
        call_cnfd.push_back(amCallCNFD(i, s_low,s_high,r, sigma, T, delta_t, delta_s1, k));
    }
    vector<vector<double>> amCallPrices{call_efd, call_ifd, call_cnfd};
    toCSV(amCallPrices, "calls.csv");
    vector<double> put_efd;
    vector<double> put_ifd;
    vector<double> put_cnfd;
    for(double i = s_low; i <= s_high; i+= delta_s1){
        put_efd.push_back(amPutEFD(i, s_low,s_high,r, sigma, T, delta_t, delta_s1, k));
        put_ifd.push_back(amPutIFD(i, s_low,s_high,r, sigma, T, delta_t, delta_s1, k));
        put_cnfd.push_back(amPutCNFD(i, s_low,s_high,r, sigma, T, delta_t, delta_s1, k));
    }
    
    vector<vector<double>> amPutPrices{put_efd, put_ifd, put_cnfd};
    toCSV(amPutPrices, "puts.csv");
    
    cout << "Graphs are in PDF." << endl;
    
    
    return 0;
}
