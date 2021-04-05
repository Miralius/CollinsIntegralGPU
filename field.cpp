#include "field.h"
extern vector<vector<complex<double>>> calculateCollinsCUDA(vector<vector<complex<double>>>& inputFunction, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v, int n1, int n2, double waveNumber, vector<double> limits, vector<vector<double>> matrixABCD);

vector<double> field::calcPoints(double begin, double end, double count) {
	auto pointValue = begin;
	vector<double> points;
	auto h = (end - begin) / count;
	for (auto i = 0; i < count; i++) {
		points.push_back(pointValue);
		pointValue += h;
	}
	return points;
}

void field::createField(const function<complex<double>(double, double)>& mode) {
	for (auto row : y) {
		calculatedField.push_back(vector<complex<double>>());
		for (auto column : x) {
			calculatedField.back().push_back(mode(column, row));
		}
	}
}

double field::maximum(vector<vector<double>>& field) {
	auto maxValue = DBL_MIN;
	for (auto row : field) {
		for (auto value : row) {
			maxValue = max(maxValue, value);
		}
	}
	return maxValue;
}

field::field() : x(vector<double>()), y(vector<double>()), calculatedField(vector<vector<complex<double>>>()) {
}

field::field(double a, double b, double n) : calculatedField(vector<vector<complex<double>>>()) {
	x = calcPoints(-a, a, n);
	y = calcPoints(-b, b, n);
	reverse(y.begin(), y.end());
}

field::field(vector<vector<complex<double>>> field, vector<double> x, vector<double> y) : calculatedField(field), x(x), y(y) {
}

vector<double> field::getX()
{
	return x;
}

vector<double> field::getY()
{
	return y;
}

vector<vector<complex<double>>> field::getCalculatedField() {
	return calculatedField;
}

void field::transform(double a, double b, double n, double wavelength, transformType transform, double z, double f) {
	auto u = calcPoints(-a, a, n);
	auto v = calcPoints(-b, b, n);
	reverse(v.begin(), v.end());
	auto k = 2 * M_PI / wavelength;
	auto limits = vector<double>({ -x.at(0), y.at(0), a, b});
	auto matrixABCD = vector<vector<double>>({ {0., 0.}, {0., 0.} });
	switch (transform) {
	case transformType::fractionalFourier:
		matrixABCD = { {cos(M_PI * z / (2 * f)), f * sin(M_PI * z / (2 * f))}, {-sin(M_PI * z / (2 * f)) / f, cos(M_PI * z / (2 * f))} };
		break;
	case transformType::fresnel:
		matrixABCD = { {1, z}, {0, 1} };
		break;
	default:
		matrixABCD = { {0., 0.}, {0., 0.} };
	}
	if (std::abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) {
		vector<vector<complex<double>>> output;
		for (auto p = 0; p < v.size(); p++) {
			output.push_back(vector<complex<double>>());
			for (auto q = 0; q < u.size(); q++) {
				output.back().push_back(matrixABCD.at(1).at(1) * calculatedField.at(p).at(q) * exp(complex<double>(0, (k * matrixABCD.at(1).at(0) * matrixABCD.at(1).at(1) * (u.at(q) * u.at(q) + v.at(p) * v.at(p))) / 2)));
			}
		}
		calculatedField = output;
	}
	else {
		calculatedField = calculateCollinsCUDA(calculatedField, x, y, u, v, x.size(), n, k, limits, matrixABCD);
	}
	x = u;
	y = v;
}


vector<vector<double>> field::abs(field& field) {
	vector<vector<double>> absField;
	for (auto row : field.getCalculatedField()) {
		absField.push_back(vector<double>());
		for (auto value : row) {
			absField.back().push_back(std::abs(value));
		}
	}
	return absField;
}

vector<vector<double>> field::arg(field& field) {
	vector<vector<double>> argField;
	for (auto row : field.getCalculatedField()) {
		argField.push_back(vector<double>());
		for (auto value : row) {
			argField.back().push_back((value.imag() < 0) ? (std::arg(value) + 2 * M_PI) : std::arg(value));
		}
	}
	return argField;
}

BMP field::createBMP(string schemeName, bool phase) {
	double minValue;
	double maxValue;
	vector<vector<double>> field;
	if (phase) {
		field = arg(*this);
		minValue = 0;
		maxValue = 2 * M_PI;
	}
	else {
		field = abs(*this);
		minValue = 0;
		maxValue = maximum(field);
	}
	vector<vector<byte>> scheme = scheme::scheme(schemeName);
	vector<vector<vector<byte>>> pixels;
	for (auto row : field) {
		pixels.push_back(vector<vector<byte>>());
		for (auto value : row) {
			pixels.back().push_back(scheme.at(static_cast<byte>(round((value - minValue) * 255 / (maxValue - minValue)))));
		}
	}
	return BMP(pixels);
}

void field::gaussMode(double sigma, double m) {
	fieldParameters = vector<double>({ sigma, m });
	auto gauss = [&](double x, double y) {
		return exp(complex<double>(-(x * x + y * y) / (fieldParameters.at(0) * fieldParameters.at(0)), 0));
	};
	createField(gauss);
}

void field::airyMode(double alpha, double beta, double alpha0, double beta0, double sigma) {
	fieldParameters = vector<double>({ alpha, beta, alpha0, beta0, sigma });
	auto airy = [&](double x, double y) {
		return exp(complex<double>(0, fieldParameters.at(0) * pow(x, 3) + fieldParameters.at(1) * pow(y, 3) + fieldParameters.at(2) * x + fieldParameters.at(3) * y));
	};
	createField(airy);
	this->transform(-x.at(0), y.at(0), x.size(), 650. / 1000000, transformType::fractionalFourier, 1000, 1000);
}

field operator*(field& left, field& right) {
	if ((left.getCalculatedField().at(1).size() != right.getCalculatedField().at(1).size()) || (left.getCalculatedField().size() != right.getCalculatedField().size())) {
		return field();
	}
	vector<vector<complex<double>>> product;
	for (auto i = 0; i < left.getCalculatedField().at(1).size(); i++) {
		product.push_back(vector<complex<double>>());
		for (auto j = 0; j < left.getCalculatedField().size(); j++) {
			product.at(i).push_back(left.getCalculatedField().at(i).at(j) * right.getCalculatedField().at(i).at(j));
		}
	}
	return field(product, left.getX(), left.getY());
}

field operator*=(field& left, field& right) {
	left = left * right;
	return left;
}
