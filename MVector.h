#ifndef MVector_h
#define MVector_h

#include <vector>

// Class that represents a mathematical vector
class MVector{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}
	MVector(std::initializer_list<double> l) : v(l) {}
	
	// access element (lvalue)
	double &operator[](int index)
	{
		if(index+1 > v.size() || index < 0){
			std::cout << "Error: Index outside range" << std::endl;
			exit(-1);
		}
		else{
			return v[index];
		}
	}
	
	// access element (rvalue)
	double operator[](int index) const {
		
		if(index+1 > v.size() || index < 0){
			std::cout << "Error: Index outside range" << std::endl;
			exit(-1);
		}
		else{
			return v[index];
		}
	}
	
	int size() const { return v.size(); } // number of elements
	
	void pushback(double x){v.push_back(x);}
	
	//Calculates the mean of a given vector of datapoints
	double average(){
		double VectorSum = 0.0;
		for(int i = 0; i < v.size(); i++){
			VectorSum += v[i];
		}
		return VectorSum / v.size();
	}
	
	
	//Calculates the moving average of a given vector, using the previous m and sucessive p data points and returns the average.
	MVector RunningAverage(int m, int p){
		MVector RunningAverageVector(v.size()-m-p);
		
		for(int i = m; i < v.size()-p; i++){
			double RunningSum = 0.0;
			for(int j = -m; j <= p; j++){
				RunningSum += v[i+j];
			}
			RunningAverageVector[i-m] = RunningSum / (m+p+1);
		}
		return RunningAverageVector;
	}
	
	//Calculates the exponential average of a given vector of datapoints with smoothing factor alpha between 0 and 1.
	MVector ExponentialAverage(double alpha){
		if(alpha < 0 || alpha > 1){
			std::cout << "ERROR: Value of alpha must be between 0 and 1 \n";
			exit(-1);
		}
		
		MVector ExponentialAverageVector(v.size());
		ExponentialAverageVector[0] = v[0];
		
		for(int i = 1; i < v.size(); i++){
			ExponentialAverageVector[i] = alpha * v[i] + (1-alpha) * ExponentialAverageVector[i-1];
		}
		return ExponentialAverageVector;
	}
	
	
	
	//Calculates the variance of the SDE modelling a stock price S following a log normal random walk
	double DailyVolatility(){
		MVector LogDataVector(v.size()-1);
		double SumOfSquares = 0.0;
		double SquareOfSum = 0.0;
		double ObservedDailyVolatility = 0.0;

		//Calculate dLnS, difference of logarithms of data values
		for(int i = 0; i < LogDataVector.size(); i++){
			LogDataVector[i] = log(v[i+1]) - log(v[i]);
		}
		
		//Calculate sum of squares of dLnS
		for(int i = 0; i < LogDataVector.size(); i++){
			SumOfSquares += pow(LogDataVector[i],2);
		}
		
		//Calculate the square of the sum of data of dLnS
		for(int i = 0; i < LogDataVector.size(); i++){
			SquareOfSum += LogDataVector[i];
		}
		SquareOfSum = pow(SquareOfSum,2);
		
		ObservedDailyVolatility = SumOfSquares/(LogDataVector.size()-1) - SquareOfSum/(LogDataVector.size()*(LogDataVector.size()-1));
		
		return ObservedDailyVolatility;
		
	}
	
private:
	std::vector<double> v;
	
};


// Overload the << operator to output MVectors to screen or file
std::ostream& operator<<(std::ostream& os, const MVector& v){
	int n = v.size();
	os << "(";
	for(int i = 0; i < n-1; i++){
		os << v[i] << ", ";
	}
	os << v[n-1] << ")" << std::endl;
	return os;
}



#endif /* MVector_h */
