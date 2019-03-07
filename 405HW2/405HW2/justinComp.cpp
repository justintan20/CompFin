//
//  justinComp.cpp
//  405HW2
//
//  Created by Justin Tan on 1/19/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justinComp.h"


using namespace std;

/*
 generates numbers from uniform distribution
 inputs: parameters a, b, m, number of numbers to generate, seed
 output: vector of generated random numbers
 */
vector<double> randUniformGen(int num, int seed){
    int a = pow(7,5);
    int m = pow(2,31) - 1;
    int b = 0;
    vector<double> result;
    long int xNow = seed;
    long int xNew = 0;
    for(int i = 0; i < num; i++){
        xNew = (a*xNow + b)%m;
        result.push_back(xNew / (double)m);
        xNow = xNew;
    }
    return result;
}

/*
 generates normal distributed random numbers using Box-Muller method
 input: vector of uniform[0,1] random numbers
 output: vector of random normal distributed numbers
 */
vector<double> normGenBM(vector<double> uniform){
    vector<double> result;
    int i = 0;
    while(i < uniform.size()){
        result.push_back(sqrt(-2*log(uniform[i]))*cos(2*M_PI*uniform[i+1]));
        result.push_back(sqrt(-2*log(uniform[i]))*sin(2*M_PI*uniform[i+1]));
        i += 2;
    }
    return result;
}

//returns mean of vector of doubles
double calcMean(vector<double> v){
    double size = v.size();
    double sum = 0.0;
    for(double i : v){
        sum += i;
    }
    return (sum / size);
}

//returns standard deviation of vector of doubles
double calcStdDev(vector<double> v){
    double size = v.size();
    double mean = calcMean(v);
    double sum = 0.0;
    for(double i : v){
        sum += pow(i-mean,2);
    }
    double stdDev = sqrt(sum / size);
    return stdDev;
}
