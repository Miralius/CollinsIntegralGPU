#pragma once
#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Windows.h>
#include <vector>
#include <complex>
#include "cuda_runtime.h"
#include "BMP.h"

constexpr auto PI = 3.1415926535897932384626433832795;

using namespace std;

inline void error(const string& s) {
	throw runtime_error(s);
}

inline void processing(int NOW, int MAX, int seconds, int timeLeft) {
	float proc, nowf, maxf;
	if (NOW == MAX) proc = 100.;
	else {
		nowf = (float)NOW;
		maxf = (float)MAX;
		proc = trunc(10000 * (nowf / maxf)) / 100;
	}
	cout << '\r' << "��������� " << setw(6) << proc << "%, ������ " << setw(6) << seconds << " ������, �������� " << setw(6) << timeLeft << " ������";
}

inline void wrongInput() {
	cout << "�������� ����! ������� ��� ���!" << endl;
	cin.clear();
	cin.ignore(cin.rdbuf()->in_avail(), '\n');
}

template <typename T> T loadingFile(string nameFile) {
	ifstream in(nameFile, ios::binary | ios::in);
	in.unsetf(ios_base::skipws);
	T data;
	if (!in.fail()) {
		if (!(in >> data)) error("���� " + nameFile + " ���� ��� �������� �������� ������!");
	}
	else {
		error("���� " + nameFile + " �� ������!");
	}
	return data;
}

template <typename T> void writingFile(T data, string nameFile) {
	ofstream output(nameFile, ios::binary | ios::trunc | ios::out);
	if (!output) {
		error("������ � ���� " + nameFile + " ����������!");
	}
	output << data;
}

vector<vector<double>> functionGauss(vector<double> x, vector<double> y, double sigma);
vector<vector<double>> functionGaussLaguerre(vector<double> x, vector<double> y, double sigma, int n, double m);

vector<vector<complex<double>>> vortex(vector<vector<double>> func, vector<double> x, vector<double> y, double n);
vector<vector<double>> abs(vector<vector<complex<double>>> field);
vector<vector<double>> arg(vector<vector<complex<double>>> field);
double minimum(vector<vector<double>> field);
double maximum(vector<vector<double>> field);

vector<vector<complex<double>>> collins(vector<vector<complex<double>>> functionVortex, vector<double> u, vector<double> v, vector<vector<double>> matrixABCD, double wavelength);
vector<vector<complex<double>>> collins(vector<vector<complex<double>>> functionVortex, vector<double> x, vector<double> y, vector<double> u, vector<double> v, vector<vector<double>> matrixABCD, double wavelength, double hx, double hy);
vector<vector<complex<double>>> cuCollins(vector<vector<complex<double>>> functionVortex, vector<double> x, vector<double> y, vector<double> u, vector<double> v, vector<vector<double>> matrixABCD, double wavelength, double hx, double hy);

#endif