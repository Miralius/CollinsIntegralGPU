#include "CollinsIntegralGPU.h"

vector<vector<complex<double>>> vortex(vector<vector<double>>& func, vector<double>& x, vector<double>& y, double n) {
	vector<vector<complex<double>>> functionVortex;
	for (int i = 0; i < y.size(); i++) {
		functionVortex.push_back(vector<complex<double>>());
		for (int j = 0; j < x.size(); j++) {
			functionVortex.at(i).push_back(func.at(i).at(j) * exp(complex<double>(0, (n * ((i < (y.size() / 2)) ? atan2(-y.at(i), x.at(j)) : (atan2(-y.at(i), x.at(j)) + 2 * PI))))));
		}
	}
	return functionVortex;
}

vector<vector<double>> abs(vector<vector<complex<double>>>& field) {
	vector<vector<double>> absField;
	for (vector<complex<double>> row : field) {
		absField.push_back(vector<double>());
		for (complex<double> value : row) {
			absField.back().push_back(abs(value));
		}
	}
	return absField;
}

vector<vector<double>> arg(vector<vector<complex<double>>>& field) {
	vector<vector<double>> argField;
	for (vector<complex<double>> row : field) {
		argField.push_back(vector<double>());
		for (complex<double> value : row) {
			argField.back().push_back((value.imag() < 0) ? (arg(value) + 2 * PI) : arg(value));
		}
	}
	return argField;
}

double minimum(vector<vector<double>>& field) {
	double minValue = DBL_MAX;
	for (vector<double> row : field) {
		for (double value : row) {
			minValue = min(minValue, value);
		}
	}
	return minValue;
}

double maximum(vector<vector<double>>& field) {
	double maxValue = DBL_MIN;
	for (vector<double> row : field) {
		for (double value : row) {
			maxValue = max(maxValue, value);
		}
	}
	return maxValue;
}