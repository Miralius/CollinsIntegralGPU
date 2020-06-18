#include "CollinsIntegralGPU.h"

int factorial(int i) {
	if (i == 0) return 1;
	else return i * factorial(i - 1);
}

double hermite(double x, int n) {
	double series = 0;
	for (int m = 0; m <= n / 2; m++) {
		series += pow(-1, m) * pow(2 * x, n - 2 * m) / (static_cast<short>(factorial(m)) * static_cast<short>(factorial(n - 2 * m)));
	}
	return series;
}

vector<vector<double>> functionGaussHermite(vector<double>& x, vector<double>& y, double sigma, int n, int m) {
	vector<vector<double>> input;
	for (int i = 0; i < y.size(); i++) {
		input.push_back(vector<double>());
		for (int j = 0; j < x.size(); j++) {
			input.at(i).push_back((hermite(y.at(i) / sigma, n) / sqrt(pow(2, n) * factorial(n) * sqrt(PI))) * (hermite(x.at(j) / sigma, m) / sqrt(pow(2, m) * factorial(m) * sqrt(PI))) * (exp(-(x.at(j) * x.at(j) + y.at(i) * y.at(i)) / (2 * sigma * sigma))));
		}
	}
	return input;
}

vector<vector<complex<double>>> superposition(vector<vector<double>> func1, vector<vector<double>> func2) {
	vector<vector<complex<double>>> summa;
	for (int i = 0; i < func1.at(1).size(); i++) {
		summa.push_back(vector<complex<double>>());
		for (int j = 0; j < func1.size(); j++) {
			summa.at(i).push_back(func1.at(i).at(j) + func2.at(i).at(j));
		}
	}
	return summa;
}

vector<vector<double>> functionGauss(vector<double>& x, vector<double>& y, double sigma) {
	vector<vector<double>> input;
	for (int i = 0; i < y.size(); i++) {
		input.push_back(vector<double>());
		for (int j = 0; j < x.size(); j++) {
			input.at(i).push_back((exp(-(x.at(j) * x.at(j) + y.at(i) * y.at(i)) / (2 * sigma * sigma))));
		}
	}
	return input;
}

vector<vector<double>> functionGaussLaguerre(vector<double>& x, vector<double>& y, double sigma, int n, double m) {
	vector<vector<double>> input;
	for (int i = 0; i < y.size(); i++) {
		input.push_back(vector<double>());
		for (int j = 0; j < x.size(); j++) {
			input.at(i).push_back((exp(-(x.at(j) * x.at(j) + y.at(i) * y.at(i)) / (2 * sigma * sigma))) * pow(((sqrt(x.at(j) * x.at(j) + y.at(i) * y.at(i))) / sigma), (int)abs(m)));
		}
	}
	return input;
}