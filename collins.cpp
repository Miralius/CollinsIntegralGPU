#include "CollinsIntegralGPU.h"

vector<vector<complex<double>>> collins(vector<vector<complex<double>>>& functionVortex, vector<double>& u, vector<double>& v, vector<vector<double>>& matrixABCD, double wavelength) {
	double k = 2 * PI / wavelength;

	vector<vector<complex<double>>> output;
	for (int i = 0; i < v.size(); i++) {
		output.push_back(vector<complex<double>>());
		for (int j = 0; j < u.size(); j++) {
			output.at(i).push_back(sqrt(matrixABCD.at(1).at(1)) * functionVortex.at(i).at(j) * exp(complex<double>(0, (k * matrixABCD.at(1).at(0) * matrixABCD.at(1).at(1) * (u.at(j) * u.at(j) + v.at(i) + v.at(i))) / 2)));
		}
	}
	return output;
}

vector<vector<complex<double>>> collins(vector<vector<complex<double>>>& functionVortex, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v, vector<vector<double>>& matrixABCD, double wavelength, double hx, double hy) {
	double k = 2 * PI / wavelength;
	int startTime(clock() / CLOCKS_PER_SEC), endTime(clock() / CLOCKS_PER_SEC), currentTime(clock() / CLOCKS_PER_SEC), progress(0);
	vector<vector<complex<double>>> output;
	for (int p = 0; p < v.size(); p++) {
		processing(++progress, (int)v.size(), clock() / CLOCKS_PER_SEC - currentTime, (endTime - startTime) * ((int)v.size() - p));
		startTime = endTime;
		output.push_back(vector<complex<double>>());
		for (int q = 0; q < u.size(); q++) {
			complex<double> value = 0;
			for (int i = 0; i < y.size(); i++) {
				for (int j = 0; j < x.size(); j++) {
					value += functionVortex.at(i).at(j) * exp(complex<double>(0, ((k / (2 * matrixABCD.at(0).at(1))) * (matrixABCD.at(0).at(0) * (y.at(i) * y.at(i) + x.at(j) * x.at(j)) - 2 * (y.at(i) * v.at(p) + x.at(j) * u.at(q)) + matrixABCD.at(1).at(1) * (v.at(p) * v.at(p) + u.at(q) * u.at(q))))));
				}
			}
			output.at(p).push_back(complex<double>(0, -(k / (2 * PI * matrixABCD.at(0).at(1)))) * value * hx * hy);
		}
		endTime = clock() / CLOCKS_PER_SEC;
	}
	return output;
}

vector<vector<complex<double>>> cuCollins(vector<vector<complex<double>>>& functionVortex, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v, vector<vector<double>>& matrixABCD, double wavelength, double hx, double hy) {
#pragma warning(push)
#pragma warning(disable:6386)
	double k = 2 * PI / wavelength;

	complex<double>** functionVortexTemp = new complex<double> * [functionVortex.size() * functionVortex.at(0).size()];
	for (int i = 0; i < functionVortex.size(); i++) {
		functionVortexTemp[i] = new complex<double>[functionVortex.size()];
		for (int j = 0; j < functionVortex.at(i).size(); j++) {
			functionVortexTemp[i][j] = functionVortex.at(i).at(j);
		}
	}

	double* xTemp = new double[x.size()];
	for (int i = 0; i < x.size(); i++) {
		xTemp[i] = x.at(i);
	}

	double* yTemp = new double[y.size()];
	for (int i = 0; i < y.size(); i++) {
		yTemp[i] = y.at(i);
	}

	double* uTemp = new double[u.size()];
	for (int i = 0; i < u.size(); i++) {
		uTemp[i] = u.at(i);
	}

	double* vTemp = new double[v.size()];
	for (int i = 0; i < v.size(); i++) {
		vTemp[i] = v.at(i);
	}

	double* constantStorage = new double[8];
	constantStorage[0] = matrixABCD.at(0).at(0);
	constantStorage[1] = matrixABCD.at(0).at(1);
	constantStorage[2] = matrixABCD.at(1).at(0);
	constantStorage[3] = matrixABCD.at(0).at(1);
	constantStorage[4] = k;
	constantStorage[5] = wavelength;
	constantStorage[6] = hx;
	constantStorage[7] = hy;


#pragma warning(pop)
	vector<vector<complex<double>>> output;
	return output;
}