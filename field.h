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
#include "scheme.h"
#include "BMP.h"

using namespace std;

enum class patternField {
	undefinedPattern, solitone, radiallySymmetricCluster
};

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
	patternField pattern; //The selected pattern of the input field.
	crossSection selectedCrossSection; //The cross section of field.
	vector<vector<double>> fieldParameters; //The parameters of input field.
	//fieldParameters.at(0): at(0) — wavelength
	//						 at(1) — N count of beams
	//						 at(2) — cluster's radius
	//						 at(3) — is the angle between the last beam and the origin zero
	//fieldParameters.at(N): at(0) — selectedInputFunction
	//						 at(1) — beam waist
	//						 at(2) — beam waist coefficient
	//						 at(3) — initial phase variant
	//						 at(4) — m / topological charge
	//						 at(5) — n
	//						 at(6) — ksi / initial transverse velocity parameter
	//						 at(7) — eta / power ratio
	//						 at(8) — initial phase variant
	vector<vector<double>> matrixABCD; //The ABCD-matrix.
	vector<vector<complex<double>>> calculatedField; //The calculated field after transform.
	vector<vector<complex<double>>> selectPatternField(vector<double>& x, vector<double>& y);
	vector<vector<complex<double>>> patternSolitone(vector<double>& x, vector<double>& y);
	vector<vector<complex<double>>> patternRadiallySymmetricCluster(vector<double>& x, vector<double>& y);
	complex<double> selectInputMode(double x, double y, int beamNumber, int iterator);
	complex<double> gaussMode(double x, double y, int beamNumber, int iterator, double beamWaistCoefficient, double initialPhaseVariant);
	complex<double> gaussHermite(double x, double y, int beamNumber);
	complex<double> gaussLaguerre(double x, double y, int beamNumber, int iterator);
	vector<vector<complex<double>>> collins(vector<vector<complex<double>>>& inputFunction, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v);
	vector<vector<complex<double>>> collinsSingular(vector<vector<complex<double>>>& inputFunction, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v);
	vector<double> calcPoints(double interval, double count, double shift);
	static double maximum(vector<vector<double>>& field);
	
public:
	field();
	field(vector<vector<complex<double>>>& superpositionField, crossSection crossSection);
	field(vector<double> limits, int n1, int n2, crossSection crossSection, vector<vector<double>> matrixABCD, patternField selectedPattern);
	static vector<vector<double>> abs(field& field);
	static vector<vector<double>> arg(field& field);
	static double power(field& field);
	void setFieldParameters(bool different, vector<vector<double>> parameters);
	vector<vector<complex<double>>> getCalculatedField();
	crossSection getCrossSection();
	BMP createBMP(string schemeName, bool phase);
	bool isSuperposition();
	bool isCalculated();
	void calculate();
};

field operator+(field& left, field& right);
vector<vector<complex<double>>> operator+(vector<vector<complex<double>>>& left, vector<vector<complex<double>>>& right);
vector<vector<complex<double>>> operator+=(vector<vector<complex<double>>>& left, vector<vector<complex<double>>>& right);

#endif