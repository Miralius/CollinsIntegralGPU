#ifndef FIELD_H
#define FIELD_H

#pragma warning(push)
#pragma warning(disable:4005)
#define _HAS_STD_BYTE 0
#pragma warning(pop)
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iomanip>
#include <functional>
#include "scheme.h"
#include "BMP.h"

using namespace std;

enum class transformType {
	fractionalFourier, fresnel
};

class field {
private:
	vector<double> x;
	vector<double> y;
	vector<vector<complex<double>>> calculatedField;
	vector<double> fieldParameters;
	vector<double> calcPoints(double begin, double end, double count);
	void createField(const function<complex<double>(double, double)>& mode);
	static double maximum(vector<vector<double>>& field);
	static double maximum(vector<complex<double>>& column);
	static vector<vector<complex<double>>> transpose(vector<vector<complex<double>>>& matrix);
public:
	field();
	field(double a, double b, double n);
	field(vector<vector<complex<double>>> field, vector<double> x, vector<double> y);
	vector<double> getX();
	vector<double> getY();
	vector<vector<complex<double>>> getCalculatedField();
	void transform(double a, double b, double n, double wavelength, transformType transformType, double z, double f=0);
	void transform(double a, double b, double n, double u, double wavelength, transformType transformType, double z_begin, double z_end, double z_n, double f=0);
	static vector<vector<double>> abs(field& field);
	static vector<vector<double>> arg(field& field);
	void normalize();
	BMP createBMP(string schemeName, bool phase);
	void gaussMode(double sigma, double m);
	void airyMode(double alpha, double beta, double alpha0, double beta0, double sigma);
};

field operator*(field& left, field& right);
field operator*=(field& left, field& right);
#endif