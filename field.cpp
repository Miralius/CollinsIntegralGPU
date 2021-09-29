#include "field.h"

std::vector<double> field::calcPoints(double begin, double end, int count) {
	auto pointValue = begin;
	std::vector<double> points;
	points.reserve(count);
	auto h = (end - begin) / count;
	for (auto i = 0; i < count; i++) {
		points.emplace_back(pointValue);
		pointValue += h;
	}
	return points;
}

void field::createField(const std::function<std::complex<double>(double, double)>& mode) {
	calculatedField.clear();
	calculatedField.reserve(y.size());
	for (auto row : y) {
		auto vectorRow = std::vector<std::complex<double>>();
		vectorRow.reserve(x.size());
		for (auto column : x) {
			vectorRow.emplace_back(mode(column, row));
		}
		calculatedField.emplace_back(vectorRow);
	}
}

double field::maximum(const std::vector<std::vector<double>>& field) {
	auto maxValue = DBL_MIN;
	for (auto &row : field) {
		for (auto value : row) {
			maxValue = max(maxValue, value);
		}
	}
	return maxValue;
}

double field::maximum(std::vector<std::complex<double>>& column) {
	auto maxValue = DBL_MIN;
	for (auto &value : column) {
		maxValue = max(maxValue, std::abs(value));
	}
	return maxValue;
}

field::field() : x(std::vector<double>()), y(std::vector<double>()), calculatedField(std::vector<std::vector<std::complex<double>>>()) {
}

field::field(double a, double b, int n) : calculatedField(std::vector<std::vector<std::complex<double>>>()) {
	x = calcPoints(a, b, n);
	y = calcPoints(a, b, n);
	reverse(y.begin(), y.end());
}

field::field(const std::vector<std::vector<std::complex<double>>>& field, const std::vector<double>& x, const std::vector<double>& y) : calculatedField(field), x(x), y(y) {
}

std::vector<double> field::getX()
{
	return x;
}

std::vector<double> field::getY()
{
	return y;
}

std::vector<std::vector<std::complex<double>>> field::getCalculatedField() {
	return calculatedField;
}

std::vector<std::vector<std::complex<double>>> field::transpose(std::vector<std::vector<std::complex<double>>>& matrix) {
	std::vector<std::vector<std::complex<double>>> transposed;
	transposed.reserve(matrix.at(0).size());
	for (auto j = 0; j < matrix.at(0).size(); j++) {
		auto column = std::vector<std::complex<double>>();
		column.reserve(matrix.size());
		for (auto i = 0; i < matrix.size(); i++) {
			column.emplace_back(matrix.at(i).at(j));
		}
		transposed.emplace_back(column);
	}
	return transposed;
}

void field::transpose() {
	std::vector<std::vector<std::complex<double>>> transposed;
	auto &field = calculatedField;
	transposed.reserve(field.at(0).size());
	for (auto j = 0; j < field.at(0).size(); j++) {
		auto column = std::vector<std::complex<double>>();
		column.reserve(field.size());
		for (auto i = 0; i < field.size(); i++) {
			column.emplace_back(field.at(i).at(j));
		}
		transposed.emplace_back(column);
	}
	calculatedField = transposed;
}

void field::normalize() {
	auto normalized = transpose(calculatedField);
	for (auto i = 0; i < normalized.size(); i++) {
		auto max = maximum(normalized.at(i));
		for (auto j = 0; j < normalized.at(0).size(); j++) {
			normalized.at(i).at(j) /= max;
		}
	}
	calculatedField = transpose(normalized);
}

void field::shift(double x, double y) {
	auto hx = this->x.at(1) - this->x.at(0);
	auto hy = this->y.at(0) - this->y.at(1);
	auto x_shift_n = static_cast<int>(x / hx);
	auto y_shift_n = static_cast<int>(y / hy);
	std::vector<std::vector<std::complex<double>>> shifted;
	shifted.reserve(calculatedField.size());
	for (unsigned long long j = 0; j < calculatedField.size(); j++) {
		auto row = std::vector<std::complex<double>>();
		row.reserve(calculatedField.at(0).size());
		for (unsigned long long i = 0; i < calculatedField.at(0).size(); i++) {
			row.emplace_back(((((i - x_shift_n) >= 0) && ((i - x_shift_n) < calculatedField.at(0).size())) && (((j - y_shift_n) >= 0) && ((j - y_shift_n) < calculatedField.size()))) ? calculatedField.at(i - x_shift_n).at(j - y_shift_n) : 0);
		}
		shifted.emplace_back(row);
	}
	calculatedField = shifted;
}

void field::rotate(double angle)
{
	std::vector<std::vector<std::complex<double>>> rotated;
	auto heightHalf = calculatedField.size() / static_cast<double>(2);
	auto widthHalf = calculatedField.at(0).size() / static_cast<double>(2);
	rotated.reserve(calculatedField.size());
	for (auto i = 0; i < calculatedField.size(); i++) {
		auto row = std::vector<std::complex<double>>();
		row.reserve(calculatedField.at(0).size());
		for (auto j = 0; j < calculatedField.at(0).size(); j++) {
			auto a = static_cast<int>(std::abs(widthHalf + (i - heightHalf) * cos(angle) + (j - widthHalf) * sin(angle)));
			auto b = static_cast<int>(std::abs(heightHalf - (i - heightHalf) * sin(angle) + (j - widthHalf) * cos(angle)));
			row.emplace_back(a < calculatedField.size() && b < calculatedField.at(0).size() ? calculatedField.at(a).at(b) : 0);
		}
		rotated.emplace_back(row);
	}
	calculatedField = rotated;
}

void field::collinsSingular(const std::vector<double>& u, const std::vector<double>& v, double wavelength, const std::vector<double>& matrixCD) {
	if (std::abs(1 + matrixCD.at(1)) < FLT_EPSILON)
	{
		this->rotate(M_PI);
	}
	std::vector<std::vector<std::complex<double>>> output;
	auto k = 2 * M_PI / wavelength;
	output.reserve(v.size());
	for (auto p = 0; p < v.size(); p++) {
		auto row = std::vector<std::complex<double>>();
		row.reserve(u.size());
		for (auto q = 0; q < u.size(); q++) {
			row.emplace_back(matrixCD.at(1) * calculatedField.at(p).at(q) * exp(std::complex<double>(0, (k * matrixCD.at(0) * matrixCD.at(1) * (u.at(q) * u.at(q) + v.at(p) * v.at(p))) / 2)));
		}
		output.emplace_back(row);
	}
	calculatedField = output;
}

bool field::hasSingularities(const std::vector<double>& z, int transformType, double parameter)
{	
	switch (transformType) {
	case 0:
		for (auto point : z)
		{
			if (std::abs(parameter * sin(M_PI * point / (2 * parameter))) < FLT_EPSILON)
			{
				return true;
			}
		}
		break;
	case 1:
		for (auto point : z)
		{
			if (std::abs(point) < FLT_EPSILON)
			{
				return true;
			}
		}
		break;
	default:
		if (std::abs(parameter) < FLT_EPSILON)
		{
			return true;
		}
	}
	return false;
}

void field::ouvFractionalFourierTransform(double a, double b, int n, double wavelength, double z, double f)
{
	auto transformType = 0;
	auto u = calcPoints(a, b, n);
	auto v = calcPoints(a, b, n);
	reverse(v.begin(), v.end());
	if (std::abs(f * sin(M_PI * z / (2 * f))) < FLT_EPSILON)
	{
		auto matrixCD = std::vector<double>({ -sin(M_PI * z / (2 * f)) / f, cos(M_PI * z / (2 * f)) });
		collinsSingular(u, v, wavelength, matrixCD);
	}
	else
	{
		auto parameters = std::vector<double>({ wavelength, z, f });
		auto dimension = std::vector<int>({ static_cast<int>(x.size()), n });
		calculatedField = calculateCollinsCUDA(calculatedField, x, y, u, v, parameters, dimension, transformType);
	}
	x = u;
	y = v;
}

void field::ouvFresnelTransform(double a, double b, int n, double wavelength, double z)
{
	auto transformType = 1;
	auto u = calcPoints(a, b, n);
	auto v = calcPoints(a, b, n);
	reverse(v.begin(), v.end());
	if (std::abs(z) < FLT_EPSILON)
	{
		auto matrixCD = std::vector<double>({ 0, 1 });
		collinsSingular(u, v, wavelength, matrixCD);
	}
	else
	{
		auto parameters = std::vector<double>({ wavelength, z });
		auto dimension = std::vector<int>({ static_cast<int>(x.size()), n });
		calculatedField = calculateCollinsCUDA(calculatedField, x, y, u, v, parameters, dimension, transformType);
	}
	x = u;
	y = v;
}

void field::ouvTransform(double a, double b, int n, double wavelength, double z, const std::vector<double>& matrixABCD)
{
	auto transformType = 2;
	auto u = calcPoints(a, b, n);
	auto v = calcPoints(a, b, n);
	reverse(v.begin(), v.end());
	if (std::abs(matrixABCD.at(1)) < FLT_EPSILON)
	{
		auto matrixCD = std::vector<double>(matrixABCD.begin() + 2, matrixABCD.end());
		collinsSingular(u, v, wavelength, matrixCD);
	}
	else
	{
		auto parameters = std::vector<double>({ wavelength, z, matrixABCD.at(0), matrixABCD.at(1), matrixABCD.at(3) });
		auto dimension = std::vector<int>({ static_cast<int>(x.size()), n });
		calculatedField = calculateCollinsCUDA(calculatedField, x, y, u, v, parameters, dimension, transformType);
	}
	x = u;
	y = v;
}

void field::ovzFractionalFourierTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u, double f)
{
	auto transformType = 0;
	auto z = calcPoints(a_z, b_z, n_z);
	auto v = calcPoints(a, b, n);
	reverse(v.begin(), v.end());
	if (hasSingularities(z, transformType, f))
	{
		throw std::runtime_error("Интервал z содержит точки вырождения интеграла Коллинза!");
	}
	else
	{
		auto parameters = std::vector<double>({ wavelength, u, f });
		auto dimension = std::vector<int>({ static_cast<int>(x.size()), n , n_z });
		calculatedField = calculateCollinsCUDA(calculatedField, x, y, z, v, parameters, dimension, transformType);
	}
	x = z;
	y = v;
}

void field::ovzFresnelTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u)
{
	auto transformType = 1;
	auto z = calcPoints(a_z, b_z, n_z);
	auto v = calcPoints(a, b, n);
	reverse(v.begin(), v.end());
	if (hasSingularities(z, transformType))
	{
		throw std::runtime_error("Интервал z содержит точки вырождения интеграла Коллинза!");
	}
	else
	{
		auto parameters = std::vector<double>({ wavelength, u });
		auto dimension = std::vector<int>({ static_cast<int>(x.size()), n , n_z });
		calculatedField = calculateCollinsCUDA(calculatedField, x, y, z, v, parameters, dimension, transformType);
	}
	x = z;
	y = v;
}

void field::ovzTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u, const std::vector<double>& matrixABCD)
{
	auto transformType = 2;
	auto z = calcPoints(a_z, b_z, n_z);
	auto v = calcPoints(a, b, n);
	reverse(v.begin(), v.end());
	if (hasSingularities(z, transformType, matrixABCD.at(1)))
	{
		throw std::runtime_error("Интервал z содержит точки вырождения интеграла Коллинза!");
	}
	else
	{
		auto parameters = std::vector<double>({ wavelength, u, matrixABCD.at(0), matrixABCD.at(1), matrixABCD.at(3) });
		auto dimension = std::vector<int>({ static_cast<int>(x.size()), n , n_z });
		calculatedField = calculateCollinsCUDA(calculatedField, x, y, z, v, parameters, dimension, transformType);
	}
	x = z;
	y = v;
}

std::vector<std::vector<double>> field::abs(const field& obj) {
	std::vector<std::vector<double>> absField;
	field data = obj;
	absField.reserve(data.getY().size());
	for (auto& row : data.getCalculatedField()) {
		auto vectorRow = std::vector<double>();
		vectorRow.reserve(data.getX().size());
		for (auto& value : row) {
			vectorRow.emplace_back(std::abs(value));
		}
		absField.emplace_back(vectorRow);
	}
	return absField;
}

std::vector<std::vector<double>> field::arg(const field& obj) {
	std::vector<std::vector<double>> argField;
	field data = obj;
	argField.reserve(data.getY().size());
	for (auto& row : data.getCalculatedField()) {
		auto vectorRow = std::vector<double>();
		vectorRow.reserve(data.getX().size());
		for (auto& value : row) {
			vectorRow.emplace_back((value.imag() < 0) ? (std::arg(value) + 2 * M_PI) : std::arg(value));
		}
		argField.emplace_back(vectorRow);
	}
	return argField;
}

BMP field::createBMP(std::string schemeName, bool phase) {
	double minValue;
	double maxValue;
	std::vector<std::vector<double>> field;
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
	std::vector<std::vector<byte>> scheme = scheme::scheme(schemeName);
	std::vector<std::vector<std::vector<byte>>> pixels;
	pixels.reserve(field.size());
	for (auto &row : field) {
		auto vectorRow = std::vector<std::vector<byte>>();
		vectorRow.reserve(field.at(0).size());
		for (auto value : row) {
			vectorRow.emplace_back(scheme.at(static_cast<byte>(round((value - minValue) * 255 / (maxValue - minValue)))));
		}
		pixels.emplace_back(vectorRow);
	}
	return BMP(pixels);
}

void field::gaussMode(double sigma, double factorSigma) {
	auto gauss = [&](double x, double y) {
		return exp(std::complex<double>(-(x * x + y * y) / (factorSigma * sigma * sigma), 0));
	};
	createField(gauss);
}

void field::gaussHermiteMode(int m, int n, double sigma)
{
	auto gaussHermite = [&](double x, double y) {
		return (std::hermite(n, y / sigma) / sqrt(pow(2, n) * tgamma(n + 1) * sqrt(M_PI))) * (std::hermite(m, x / sigma) / sqrt(pow(2, m) * tgamma(m + 1) * sqrt(M_PI))) * exp(std::complex<double>(-(x * x + y * y) / (2 * sigma * sigma), 0));
	};
	createField(gaussHermite);
}

void field::gaussLaguerreMode(int m, int n, double sigma)
{
	auto gaussLaguerre = [&](double x, double y) {
		auto gaussMode = exp(std::complex<double>(-(x * x + y * y) / (sigma * sigma), 0));
		auto phi = y > 0 ? atan2(y, x) : atan2(y, x) + 2 * M_PI;
		auto vortexMode = exp(std::complex<double>(0, m * phi));
		return (1 / sigma) * sqrt(2 * tgamma(n + 1) / (M_PI * tgamma(n + std::abs(m) + 1))) * gaussMode * pow(((sqrt(2 * (x * x + y * y))) / sigma), std::abs(m)) * std::assoc_laguerre(n, std::abs(m), 2 * (x * x + y * y) / (sigma * sigma)) * vortexMode;
	};
	createField(gaussLaguerre);
}

void field::airyMode(double alpha, double beta, double alpha0, double beta0) {
	auto airy = [&](double x, double y) {
		return exp(std::complex<double>(0, alpha * pow(x, 3) + beta * pow(y, 3) + alpha0 * x + beta0 * y));
	};
	createField(airy);
	this->ouvFractionalFourierTransform(y.back(), -y.back(), static_cast<int>(y.size()), 650. / 1000000, 1000, 1000);
}

void field::vortexMode(double m)
{
	auto vortex = [&](double x, double y) {
		auto phi = y > 0 ? atan2(y, x) : atan2(y, x) + 2 * M_PI;
		return exp(std::complex<double>(0, m * phi));
	};
	createField(vortex);
}

void field::initialTransverseVelocityAndPowerFactorExpMode(double ksi, double eta, double sigma, double x_shift, double y_shift) {
	auto initialTransverseVelocityAndPowerFactorExp = [&](double x, double y) {
		return exp(std::complex<double>(0, -ksi * (sqrt(eta) / (sigma * sigma)) * (y_shift * x - x_shift * y)));
	};
	createField(initialTransverseVelocityAndPowerFactorExp);
}

field operator+(const field& left, const field& right) {
	field leftField = left;
	field rightField = right;
	auto leftMatrix = leftField.getCalculatedField();
	auto rigthMatrix = rightField.getCalculatedField();
	if ((leftMatrix.at(0).size() != rigthMatrix.at(0).size()) || (leftMatrix.size() != rigthMatrix.size())) {
		return field();
	}
	std::vector<std::vector<std::complex<double>>> sum;
	sum.reserve(leftMatrix.size());
	for (auto i = 0; i < leftMatrix.size(); i++) {
		auto row = std::vector<std::complex<double>>();
		row.reserve(leftMatrix.at(0).size());
		for (auto j = 0; j < leftMatrix.at(0).size(); j++) {
			row.emplace_back(leftMatrix.at(i).at(j) + rigthMatrix.at(i).at(j));
		}
		sum.emplace_back(row);
	}
	return field(sum, leftField.getX(), leftField.getY());
}

field operator+=(field& left, const field& right) {
	left = left + right;
	return left;
}

field operator*(const field& left, const field& right) {
	field leftField = left;
	field rightField = right;
	auto leftMatrix = leftField.getCalculatedField();
	auto rigthMatrix = rightField.getCalculatedField();
	if ((leftMatrix.at(0).size() != rigthMatrix.at(0).size()) || (leftMatrix.size() != rigthMatrix.size())) {
		return field();
	}
	std::vector<std::vector<std::complex<double>>> product;
	product.reserve(leftMatrix.size());
	for (auto i = 0; i < leftMatrix.size(); i++) {
		auto row = std::vector<std::complex<double>>();
		product.reserve(leftMatrix.at(0).size());
		for (auto j = 0; j < leftMatrix.at(0).size(); j++) {
			row.emplace_back(leftMatrix.at(i).at(j) * rigthMatrix.at(i).at(j));
		}
		product.emplace_back(row);
	}
	return field(product, leftField.getX(), leftField.getY());
}

field operator*=(field& left, const field& right) {
	left = left * right;
	return left;
}
