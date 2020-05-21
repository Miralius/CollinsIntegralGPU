#include "CollinsIntegralGPU.h"

vector<vector<double>> functionGauss(vector<double> x, vector<double> y, double sigma) {
	vector<vector<double>> input;
	for (int i = 0; i < y.size(); i++) {
		input.push_back(vector<double>());
		for (int j = 0; j < x.size(); j++) {
			input.at(i).push_back((exp(-(x.at(j) * x.at(j) + y.at(i) * y.at(i)) / (2 * sigma * sigma))));
		}
	}
	return input;
}

vector<vector<double>> functionGaussLaguerre(vector<double> x, vector<double> y, double sigma, int n, double m) {
	vector<vector<double>> input;
	for (int i = 0; i < y.size(); i++) {
		input.push_back(vector<double>());
		for (int j = 0; j < x.size(); j++) {
			input.at(i).push_back((exp(-(x.at(j) * x.at(j) + y.at(i) * y.at(i)) / (2 * sigma * sigma))) * pow(((sqrt(x.at(j) * x.at(j) + y.at(i) * y.at(i))) / sigma), (int)abs(m)));
		}
	}
	return input;
}