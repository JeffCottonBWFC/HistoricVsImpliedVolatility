#include <iostream>
#include <fstream>
#include <cmath>
#include "MVector.h"



//Payoff function for a call option using transformed variables
double payoff(double y, double X){
	return fmax(X*exp(y) - X, 0);
}

//Function A(x)
double A(double x, double r, double sigma, double T, double k){
	return exp(-0.5*k*x - 0.125*pow(sigma,2)*pow(k,2)*T - r*T)/sqrt(2*pow(sigma,2)*M_PI*T);
}

//Function B(x,y)
double B(double x, double y, double sigma, double T, double k){
	return exp(-pow((x-y),2)/(2*pow(sigma,2)*T) + 0.5*k*y);
}



// integrand function
double f(double x, double y, double X, double r, double sigma, double T, double k){
	return B(x, y, sigma, T, k) * payoff(y,X);
}


// function to implement Simpson's integration
double integrate(double a, double b, int n, double x, double X, double r, double sigma, double T, double k){
	
	//Check number of divisions is even (as required by Simpson's algorithm)
	if(n % 2 != 0){
		std::cout << "The number of divisions 'n' must be even \n";
		exit(-1);
	}
		
	//Simpsons integration algorithm
	double sum = 0.0;
	double h = (b-a)/n;
	double EvensSum = 0.0;
	double OddsSum = 0.0;

	for(int i = 1; i <= n/2 - 1; i++){
		EvensSum += f(x, a+2*i*h, X, r, sigma, T, k);
	}
	
	for(int i = 1; i <= n/2; i++){
		OddsSum += f(x, a+(2*i-1)*h, X, r, sigma, T, k);
	}
		
	sum = h/3 * (f(x, a, X, r, sigma, T, k) + 2 * EvensSum + 4 * OddsSum + f(x, b, X, r, sigma, T, k));
		
	return sum;
}




//Function g(x) for Newton-Secant root finding algorithm
double g(double x){
	return pow(x,2) - 2;
}

//Function which takes in two initial guesses, a tolerance of convergence and the maximum number of iterations and returns a boolean of convergence and the convergent value of x
bool NewtonSecant(double &result, double x0, double x1, double tol, int maxiter){

	int i = 0;
	
	for (i = 0; i < maxiter; i++){
		result = x1 - g(x1)*(x1-x0)/(g(x1)-g(x0));
		
		if(std::abs(g(result) - g(x1)) < tol){
			break;
		}
			
		x0 = x1;
		x1 = result;
	}
	
	if(i == maxiter){
		std::cout << "Max iterations reached without convergence \n";
		return false;
	}else{
		std::cout << "The convergent value of the given function is " << result << "\n";
		return true;
	}
	
}


int main() {
	

}



//OLD CODE


	//Newton Secant Example Code
//double result = 0.0;
//double x0 = 0;
//double x1 = 5;
//double tol = 1e-5;
//int maxiter = 100;
//
//NewtonSecant(result, x0, x1, tol, maxiter);

	//Evaluating price of European call options usings procedural approach
//	Test data values
//	double X = 100;
//	double r = 0.06;
//	double sigma = 0.2;
//	double T = 0.75;
//	double S_0 = 97;
//
//	double k = 2*r/pow(sigma,2) - 1;
//	double x = log(S_0/X);
//
//	std::cout << integrate(-10, 10, 1000, x, X, r, sigma, T, k) << "\n";

	//Simpsons rule integration varying T from 0 to 1 in increments of 0.25
// Test data values
//	double X = 10;
//	double r = 0.05;
//	double sigma = 0.1;
//	double T = 0;
//
//	while(T <= 1){
//		std::cout << std::setprecision(7) << integrate(0,1,6, X, r, sigma, T) << "\n";
//		T+=0.25;
//	}


	//Tesla closing prices reading from file test
//	MVector TeslaClosingPrices;
//	long double inputVal = 0.0;
//
//	//Read TeslaClosingPrices from file
//	std::ifstream input("TeslaClosingPrices.txt");
//	if (!input){
//		std::cout << "Could not open file for reading" << std::endl;
//	}
//	else{
//		while(!input.eof()){
//			input >> inputVal;
//			TeslaClosingPrices.pushback(inputVal);
//		}
//	}


