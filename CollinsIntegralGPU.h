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

enum class scheme {
	black_white, red, fire
};

inline void error(const string& s) {
	throw runtime_error(s);
}

inline void processing(int now, int max, int seconds, int timeLeft) {
	double percent;
	if (now == max) percent = static_cast<double>(100);
	else {
		percent = trunc(10000 * (static_cast<double>(now) / static_cast<double>(max))) / 100;
	}
	cout << '\r' << "Выполнено " << setw(6) << percent << "%, прошло " << setw(6) << seconds << " секунд, осталось " << setw(6) << timeLeft << " секунд";
}

inline void wrongInput() {
	cout << "Неверный ввод! Введите ещё раз!" << endl;
	cin.clear();
	cin.ignore(cin.rdbuf()->in_avail(), '\n');
}

template <typename T> T loadingFile(string nameFile) {
	ifstream in(nameFile, ios::binary | ios::in);
	in.unsetf(ios_base::skipws);
	T data;
	if (!in.fail()) {
		if (!(in >> data)) error("Файл " + nameFile + " пуст или содержит неверные данные!");
	}
	else {
		error("Файл " + nameFile + " не найден!");
	}
	return data;
}

template <typename T> vector<T> loadingData(string nameFile) {
	ifstream in(nameFile);
	vector<T> vectorName;
	if (!in.fail()) {
		T buffer;
		while (in >> buffer) {
			if (in.eof()) break;
			vectorName.push_back(buffer);
		}
		in.close();
		if (vectorName.size() == 0) error("Файл " + nameFile + " пуст или содержит неверные данные!");
	}
	else error("Файл " + nameFile + " не найден!");
	return vectorName;
}

template <typename T> void writingFile(T data, string nameFile) {
	ofstream output(nameFile, ios::binary | ios::trunc | ios::out);
	if (!output) {
		error("Запись в файл " + nameFile + " невозможна!");
	}
	output << data;
}

vector<vector<double>> functionGauss(vector<double>&  x, vector<double>&  y, double sigma);
vector<vector<double>> functionGaussLaguerre(vector<double>&  x, vector<double>&  y, double sigma, int n, double m);

vector<vector<complex<double>>> vortex(vector<vector<double>>&  func, vector<double>&  x, vector<double>&  y, double n);
vector<vector<double>> abs(vector<vector<complex<double>>>&  field);
vector<vector<double>> arg(vector<vector<complex<double>>>&  field);
double minimum(vector<vector<double>>&  field);
double maximum(vector<vector<double>>&  field);
vector<vector<unsigned char>> applyScheme(scheme schemeName);
vector<vector<vector<unsigned char>>> fieldToBMP(vector<vector<double>> field, scheme schemeName);

vector<vector<complex<double>>> collins(vector<vector<complex<double>>>&  functionVortex, vector<double>& u, vector<double>&  v, vector<vector<double>>& matrixABCD, double wavelength);
vector<vector<complex<double>>> collins(vector<vector<complex<double>>>& functionVortex, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v, vector<vector<double>>& matrixABCD, double wavelength, double hx, double hy);
vector<vector<complex<double>>> cuCollins(vector<vector<complex<double>>>& functionVortex, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v, vector<vector<double>>& matrixABCD, double wavelength, double hx, double hy);

#endif