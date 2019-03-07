//
//  main.cpp
//  test
//
//  Created by Justin Tan on 1/10/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <fstream>
#include <chrono>

using namespace std;

/*
 generates numbers from uniform distribution
 inputs: parameters a, b, m, number of numbers to generate, seed
 output: vector of generated random numbers
*/
vector<double> randUniformGen(int a, int b, int m, int num, int seed){
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
 generates random numbers from discrete distribution
 inputs: vector of uniformly distibuted random numbers, vector of range of probabilities, vector of correspoding values
 output: vector of generated random discrete numbers
*/
vector<double> randDiscreteGen(vector<double> uniform, vector<double> pRange, vector<double> x){
    vector<double> result;
    for(int i = 0; i < uniform.size(); i++){
        for(int j = 0; j < pRange.size(); j++){
            if(uniform[i] <= pRange[j]){
                result.push_back(x[j]);
                break;
            }
        }
    }
    return result;
}

/*
 generates random binomial distributed numbers
 inputs: parameters n and p, number of numbers to generate
 output: vector of generated numbers
 */
vector<double> binomGen(int n, double p, int num){
    vector<double> result;
    int a = pow(7,5);
    int m = pow(2,31) - 1;
    int b = 0;
    int seed = 300;
    vector<double> uniform = randUniformGen(a, b, m, n*num, seed);
    vector<int> bernoulli;
    for(double i : uniform){
        if(i <= p){
            bernoulli.push_back(1);
        }
        else{
            bernoulli.push_back(0);
        }
    }
    int index = -1;
    for(int j = 0; j < bernoulli.size(); j++){
        if(j % 44 == 0){
            result.push_back(bernoulli[j]);
            index++;
        }
        else{
            result[index] += bernoulli[j];
        }
    }
    return result;
}

/*
 generates random exponential distribution numbers
 inputs: vector of uniformly distributed numbers, lambda
 output: vector of exponentialy distributed numbers
 */
vector<double> expGen(vector<double> uniform, double lambda){
    vector<double> result;
    for(double i : uniform){
        result.push_back((-1.0)*lambda*log(1-i));
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

/*
 generates normal distributed random numbers using Polar-Marsaglia method
 input: vector of uniform[0,1] random numbers
 output: vector of random normal distributed numbers
 */
vector<double> normGenPM(vector<double> uniform){
    vector<double> result;
    int i = 0;
    while(i < uniform.size()){
        double v1 = 2*uniform[i] - 1;
        double v2 = 2*uniform[i+1] - 1;
        double w = v1*v1 + v2*v2;
        if(w <= 1){
            result.push_back(sqrt(((-2)*log(w))/w)*v1);
            result.push_back(sqrt(((-2)*log(w))/w)*v2);
        }
        i += 2;
    }
    return result;
}

vector<double> normGenPM2(vector<double> uniform, int n){
    vector<double> result;
    int i = 0;
    while(result.size() < n){
        double v1 = 2*uniform[i] - 1;
        double v2 = 2*uniform[i+1] - 1;
        double w = v1*v1 + v2*v2;
        if(w <= 1){
            result.push_back(sqrt(((-2)*log(w))/w)*v1);
            result.push_back(sqrt(((-2)*log(w))/w)*v2);
        }
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

int main(int argc, const char * argv[]) {
    //uniform setup
    int a = pow(7,5);
    int m = pow(2,31) - 1;
    int b = 0;
    int seed = 80000;
    int n = 10000;
    vector<double> nums;
    cout << "Generating 10,000 U[0,1] random numbers...";
    nums = randUniformGen(a, b, m, n, seed);
    cout << "Done!" << endl;
    cout << "Mean: " << calcMean(nums) << "\tStd. Dev.: " << calcStdDev(nums) << endl;
    //built-in uniform
    cout << "Generating 10,000 U[0,1] random numbers using built-in function...";
    default_random_engine generator;
    uniform_real_distribution<double> dist(0.0,1.0);
    vector<double> builtIn;
    for(int i = 0; i < 10000; i++){
        builtIn.push_back(dist(generator));
    }
    cout << "Done!" << endl;
    cout << "Mean: " << calcMean(builtIn) << "\tStd. Dev.: " << calcStdDev(builtIn) << endl;
    cout << "Comparing the 2, we see that the built-in random numbers have marginally larger mean and standard deviation, but barely noticable.\n" << endl;
    cout << "Generating 10,000 discrete distribution random numbers...";
    //discrete
    vector<double> x {-1, 0, 1, 2};
    vector<double> pRanges {0.3, 0.3+0.35, 0.3+0.35+0.2, 0.3+0.35+0.2+0.15};
    vector<double> discreteRand = randDiscreteGen(nums, pRanges, x);
    cout << "Done!" << endl;
    ofstream disHist("discrete.txt");
    for(int i = -1; i < 3; i++){
        disHist << i << ":\t";
        int count = 0;
        for(double j : discreteRand){
            if(j == i){
                count++;
                if(count == 50){
                    disHist << "*";
                    count = 0;
                }
            }
        }
        disHist << endl;
    }
    disHist.close();
    cout << "Histogram is in directory as file: discrete.txt" << endl;
    cout << "Mean: " << calcMean(discreteRand) << "\tStd. Dev.: " << calcStdDev(discreteRand) << "\n" << endl;
    
    cout << "Generating 1,000 binomial(n=44, p=0.64) distribution random numbers...";
    //binomial
    vector<double> binomials = binomGen(44, 0.64, 1000);
    cout << "Done!" << endl;
    ofstream histogram("binomial.txt");
    for(int i = 0; i < 45; i++){
        histogram << i << ":\t";
        for(double b : binomials){
            if(b == i){
                histogram << "*";
            }
        }
        histogram << endl;
    }
    histogram.close();
    cout << "Histogram is in directory as file: binomial.txt" << endl;
    int count40 = 0;
    for(double i : binomials){
        if(i >= 40){
            count40++;
        }
    }
    double prob40 = ((double)count40)/binomials.size();
    cout << "P(X >= 40): " << prob40 << endl;
    cout << "From online resources, the exact calculated probability of P(X >= 40) is around 0.00008, which is very close to our result, 0, as well. If we simulate more than 1,000 random binomial numbers we should approach this result.\n" << endl;
    
    cout << "Generating 10,000 exponentially distributed random numbers (lambda = 1.5)...";
    //exponential
    vector<double> expNums = expGen(nums, 1.5);
    cout << "Done!" << endl;
    int count1 = 0;
    int count4 = 0;
    for(double i : expNums){
        if(i >= 1.0){
            count1++;
        }
        if(i >= 4.0){
            count4++;
        }
    }
    double size = expNums.size();
    double prob1 = count1 / size;
    double prob4 = count4 / size;
    cout << "P(X >= 1): " << prob1 << "\tP(X >= 4): " << prob4 << endl;
    cout << "Mean: " << calcMean(expNums) << "\tStd. Dev.: " << calcStdDev(expNums) << endl;
    ofstream expHist("exp.txt");
    double low = 0;
    double high = 0.1;
    double max = 0;
    for(double num : expNums){
        if(num > max){
            max = num;
        }
    }
    double intervals = max/0.1;
    for(int i = 0; i <= intervals; i++){
        expHist << low << ":\t ";
        int count = 0;
        for(double j : expNums){
            if(j > low && j <= high){
                count++;
                if(count == 5){
                    expHist << "*";
                    count = 0;
                }
            }
        }
        expHist << endl;
        low = high;
        high = low + 0.1;
    }
    expHist.close();
    cout << "Histogram is in directory as file: exp.txt\n" << endl;
    //normal
    cout << "Generating 5,000 U[0,1] random numbers...";
    vector<double> num5 = randUniformGen(a, b, m, 5000, seed);
    cout << "Done!\n" << endl;
    //box-muller
    cout << "Generating 5,000 N[0,1] random numbers (Box-Muller)...";
    auto start = chrono::high_resolution_clock::now();
    vector<double> normNumBM = normGenBM(num5);
    auto stop = chrono::high_resolution_clock::now();
    cout << "Done!" << endl;
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Mean: "<< calcMean(normNumBM) << "\tStd.Dev.: " << calcStdDev(normNumBM) << "\tExecution time (microseconds): " << duration.count() <<endl;
    //polar-marsaglia
    cout << "\nGenerating N[0,1] random numbers (Polar-Marsaglia)...";
    start = chrono::high_resolution_clock::now();
    vector<double> normNumPM = normGenPM(num5);
    stop = chrono::high_resolution_clock::now();
    cout << "Done!" << endl;
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Mean: "<< calcMean(normNumPM) << "\tStd.Dev.: " << calcStdDev(normNumPM) << "\tExecution time (microseconds): " << duration.count() <<endl;
    vector<double> timeTestUniform = randUniformGen(a, b, m, 8000, seed);
    start = chrono::high_resolution_clock::now();
    cout << "\nGenerating 5,000 N[0,1] random numbers (Polar-Marsaglia)...";
    vector<double> timeTestPM = normGenPM2(timeTestUniform, 5000);
    stop = chrono::high_resolution_clock::now();
    cout << "Done!" << endl;
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Mean: "<< calcMean(timeTestPM) << "\tStd.Dev.: " << calcStdDev(timeTestPM) << "\tExecution time (microseconds): " << duration.count() <<endl;
        cout << "\nBy comparing the execution times, we see that the Polar-Marsaglia method is more efficient.\n" << endl;
    return 0;
}
