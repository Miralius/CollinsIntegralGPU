﻿#ifndef FIELD_H
#define FIELD_H

#define _USE_MATH_DEFINES
#define _HAS_STD_BYTE 0
#include <cmath>
#include <complex>
#include <iomanip>
#include "scheme.h"
#include "BMP.h"

using namespace std;

enum class inputField {
	undefined, gauss, gaussHermite, gaussLaguerre
};

enum class crossSection {
	Oxy, Oxz, Oyx, Oyz, Ozx, Ozy
};

class field {
private:
	int n1; //Dimension of the input field.
	int n2; //Dimension of the output field.
	bool superposition; //The value of superposition attribute.
	bool calculated; //This value means field is calculated.
	vector<double> limits; //The limits of the integration field.
	inputField inputFunction; //The type of input function.
	crossSection selectedCrossSection; //The cross section of field.
	vector<double> fieldParameters; //The parameters of input field.
	vector<vector<double>> matrixABCD; //The ABCD-matrix.
	vector<vector<complex<double>>> calculatedField; //The calculated field after transform.
	vector<vector<complex<double>>> gauss(vector<double>& x1, vector<double>& x2, double sigma, double m);
	vector<vector<complex<double>>> gaussHermite(vector<double>& x1, vector<double>& x2, double sigma, double m, double n);
	vector<vector<complex<double>>> gaussLaguerre(vector<double>& x1, vector<double>& x2, double sigma, double m, double n);
	vector<vector<complex<double>>> selectInputField(vector<double>& x1, vector<double>& x2);
	vector<vector<complex<double>>> collins(vector<vector<complex<double>>>& inputFunction, vector<double>& x1, vector<double>& x2, vector<double>& x3, vector<double>& x4);
	vector<vector<complex<double>>> collinsSingular(vector<vector<complex<double>>>& inputFunction, vector<double>& x1, vector<double>& x2, vector<double>& x3, vector<double>& x4);
	vector<double> calcPoints(double interval, double count);
	static double maximum(vector<vector<double>>& field);
	void processing(int now, int max, int seconds, int timeLeft);
	
public:
	field();
	field(vector<vector<complex<double>>>& superpositionField, crossSection crossSection);
	field(vector<double> limits, int n1, int n2, crossSection crossSection, vector<vector<double>> matrixABCD, inputField inputFunction, vector<double> fieldParameters);
	static vector<vector<double>> abs(field& field);
	static vector<vector<double>> arg(field& field);
	vector<vector<complex<double>>> getCalculatedField();
	crossSection getCrossSection();
	BMP createBMP(string schemeName, bool phase);
	bool isSuperposition();
	bool isCalculated();
	void calculate();
};

field operator+(field& firstField, field& secondField);

#endif