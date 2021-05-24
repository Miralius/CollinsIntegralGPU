#ifndef FIELD_H
#define FIELD_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iomanip>
#include <functional>
#include "scheme.h"
#include "BMP.h"

//Функция вычисления интеграла Коллинза с помощью технологии CUDA. Входные параметры: input — входное поле, x1 — вектор x, x2 — вектор y, x3 — вектор u или z, x4 — вектор v, parameters — массив параметров преобразования, dimension — массив размерностей, transformType — вид преобразования: 0 — дробное преобразование Фурье, 1 — преобразование Френеля, 2 — другое (с заданной определенной ABCD матрицей). Выходной параметр: output — выходное поле (результат вычисления интеграла Коллинза).
extern std::vector<std::vector<std::complex<double>>> calculateCollinsCUDA(const std::vector<std::vector<std::complex<double>>>& input, const std::vector<double>& x1, const std::vector<double>& x2, const std::vector<double>& x3, const std::vector<double>& x4, const std::vector<double>& parameters, const std::vector<int>& dimension, int transformType);

//Класс светового поля
class field {
private:
	std::vector<double> x; // Вектор x (до преобразования) или u (после преобразования), или z (в случае продольного сечения преобразования)
	std::vector<double> y; // Вектор y (до преобразования) или v (после преобразования)
	std::vector<std::vector<std::complex<double>>> calculatedField;
	static std::vector<double> calcPoints(double begin, double end, int count);
	void createField(const std::function<std::complex<double>(double, double)>& mode);
	static double maximum(const std::vector<std::vector<double>>& field);
public:
	field();
	field(double a, double b, int n);
	field(const std::vector<std::vector<std::complex<double>>>& field, const std::vector<double>& x, const std::vector<double>& y);
	std::vector<double> getX();
	std::vector<double> getY();
	std::vector<std::vector<std::complex<double>>> getCalculatedField();
	void shift(double x, double y);
	void rotate(double angle);
	void collinsSingular(const std::vector<double>& u, const std::vector<double>& v, double wavelength, const std::vector<double>& matrixCD);
	static bool hasSingularities(const std::vector<double>& z, int transformType, double parameter = 0);
	void ouvFractionalFourierTransform(double a, double b, int n, double wavelength, double z, double f);
	void ouvFresnelTransform(double a, double b, int n, double wavelength, double z);
	void ouvTransform(double a, double b, int n, double wavelength, double z, const std::vector<double>& matrixABCD);
	void ovzFractionalFourierTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u, double f);
	void ovzFresnelTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u);
	void ovzTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u, const std::vector<double>& matrixABCD);
	static std::vector<std::vector<double>> abs(const field& field);
	static std::vector<std::vector<double>> arg(const field& field);
	BMP createBMP(std::string schemeName, bool phase);
	void gaussMode(double sigma, double factorSigma);
	void airyMode(double alpha, double beta, double alpha0, double beta0);
	void initialTransverseVelocityAndPowerFactorExpMode(double ksi, double eta, double sigma, double x_shift, double y_shift);
};

field operator+(const field& left, const field& right);
field operator+=(field& left, const field& right);
field operator*(const field& left, const field& right);
field operator*=(field& left, const field& right);
#endif