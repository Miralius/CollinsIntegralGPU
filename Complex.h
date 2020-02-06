#pragma once

class Complex {
private:
	double value[2];

public:
	double real() const {
		return value[0];
	}

	double imag() const {
		return value[1];
	}

	Complex() {
		value[0] = 0;
		value[1] = 0;
	}

	Complex(double obj1, double obj2) {
		value[0] = obj1;
		value[1] = obj2;
	}

	Complex& operator=(const Complex right) {
		value[0] = right.real();
		value[1] = right.imag();
		return *this;
	}

	Complex& operator=(const double right) {
		value[0] = right;
		value[1] = 0;
		return *this;
	}
};

Complex operator*(const double& left, const Complex& right) {
	return Complex(right.real() * left, right.imag() * left);
}