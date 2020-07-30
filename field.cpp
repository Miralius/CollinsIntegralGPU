#include "field.h"
extern vector<vector<complex<double>>> calculateCollinsCUDA(vector<vector<complex<double>>>& inputFunction, vector<double>& x1, vector<double>& x2, vector<double>& x3, vector<double>& x4, int n1, int n2, double waveNumber, vector<double> limits, vector<vector<double>> matrixABCD);

vector<vector<complex<double>>> field::gauss(vector<double>& x1, vector<double>& x2, double sigma, double m) {
	vector<vector<complex<double>>> input;
	for (auto i = 0; i < x2.size(); i++) {
		input.push_back(vector<complex<double>>());
		for (auto j = 0; j < x1.size(); j++) {
			auto phi = (i < (x2.size() / 2)) ? atan2(x2.at(i), -x1.at(j)) : (atan2(x2.at(i), -x1.at(j)) + 2 * M_PI);
			input.back().push_back(exp(complex<double>(-(x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (2 * sigma * sigma), 0)) * exp(complex<double>(0, m * phi)));
		}
	}
	return input;
}

vector<vector<complex<double>>> field::gaussHermite(vector<double>& x1, vector<double>& x2, double sigma, double m, double n) {
	vector<vector<complex<double>>> input;
	for (auto i = 0; i < x2.size(); i++) {
		input.push_back(vector<complex<double>>());
		for (auto j = 0; j < x1.size(); j++) {
			//std::hermite() requires C++17
			input.back().push_back((hermite(static_cast<int>(n), x2.at(i) / sigma) / sqrt(pow(2, n) * tgamma(n + 1) * sqrt(M_PI))) * (hermite(static_cast<int>(m), x1.at(j) / sigma) / sqrt(pow(2, m) * tgamma(m + 1) * sqrt(M_PI))) * (exp(-(x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (2 * sigma * sigma))));
		}
	}
	return input;
}

vector<vector<complex<double>>> field::gaussLaguerre(vector<double>& x1, vector<double>& x2, double sigma, double m, double n) {
	vector<vector<complex<double>>> input;
	for (auto i = 0; i < x2.size(); i++) {
		input.push_back(vector<complex<double>>());
		for (auto j = 0; j < x1.size(); j++) {
			//std::assoc_laguerre() requires C++17
			auto phi = (i < (x2.size() / 2)) ? atan2(x2.at(i), -x1.at(j)) : (atan2(x2.at(i), -x1.at(j)) + 2 * M_PI);
			input.back().push_back((1 / sigma) *sqrt(2 * tgamma(n + 1) / (M_PI * tgamma(n + std::abs(m) + 1))) * exp(-(x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (sigma * sigma)) * pow(((sqrt(2 * (x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)))) / sigma), std::abs(m)) * assoc_laguerre(static_cast<unsigned int>(n), static_cast<unsigned int>(std::abs(m)), 2 * (x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (sigma * sigma)) * exp(complex<double>(0, m * phi)));
		}
	}
	return input;
}

vector<vector<complex<double>>> field::gaussWang2008(vector<double>& x1, vector<double>& x2, double wavelength, double sigma, double m, double L0, double N, double k) {
	vector<vector<complex<double>>> input;
	for (auto i = 0; i < x2.size(); i++) {
		input.push_back(vector<complex<double>>());
		for (auto j = 0; j < x1.size(); j++) {
			input.back().push_back(exp(complex<double>(0, 2 * M_PI / wavelength * L0)) * exp(complex<double>(-(x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (sigma * sigma), 0)) * exp(complex<double>(0, m * M_PI * (2 * k - 1) / N)));
		}
	}
	return input;
}

vector<vector<complex<double>>> field::gaussSong2018(vector<double>& x1, vector<double>& x2, double wavelength, double sigma, double n0, double radius, double velocity, double relativePower, double beamNumber, double countOfBeams) {
	vector<vector<complex<double>>> input;
	for (auto i = 0; i < x2.size(); i++) {
		input.push_back(vector<complex<double>>());
		for (auto j = 0; j < x1.size(); j++) {
			auto phi_n = 2 * M_PI * beamNumber / countOfBeams;
			auto c_xn = radius * cos(phi_n);
			auto c_yn = radius * sin(phi_n);
			input.back().push_back(exp(complex<double>(-((x1.at(j) - c_xn) * (x1.at(j) - c_xn) + (x2.at(i) - c_yn) * (x2.at(i) - c_yn)) / (2 * sigma * sigma), 0)) / (2 * M_PI * sigma * sigma));
		}
	}
	double P0;
	if (fieldParameters.at(7) == 0) {
		calculatedField = input;
		P0 = power(*this);
		fieldParameters.at(7) = P0;
	}
	else {
		P0 = fieldParameters.at(7);
	}
	auto k = 2 * M_PI * n0 / wavelength;
	auto gamma = sqrt(relativePower) / (sqrt(P0) * k * sigma * sigma);
	auto beta = gamma * sqrt(P0);
	for (auto i = 0; i < x2.size(); i++) {
		for (auto j = 0; j < x1.size(); j++) {
			auto phi_n = 2 * M_PI * beamNumber / countOfBeams;
			input.at(i).at(j) *= sqrt(P0) / (sqrt(M_PI) * sigma) * exp(complex<double>(0, - k * velocity * beta * radius * (sin(phi_n) * x1.at(j) - cos(phi_n) * x2.at(i))));
		}
	}
	return input;
}

vector<vector<complex<double>>> field::clusterModeGaussWang2008(double wavelength, double sigma, double m, double L0, double N, double radius) {
	vector<vector<complex<double>>> inputField;
	for (auto i = 1; i <= static_cast<int>(N); i++) {
		auto x = calcPoints(limits.at(0), n1, radius * cos(M_PI * (2 * static_cast<long long>(i) - 1) / N));
		auto y = calcPoints(limits.at(1), n1, radius * sin(M_PI * (2 * static_cast<long long>(i) - 1) / N));
		reverse(y.begin(), y.end());
		inputField += gaussWang2008(x, y, wavelength, sigma, m, L0, N, static_cast<double>(i));
	}
	return inputField;
}

vector<vector<complex<double>>> field::clusterModeGaussSong2018(vector<double>& x1, vector<double>& x2, double wavelength, double sigma, double n0, double radius, double velocity, double relativePower, double countOfBeams) {
	vector<vector<complex<double>>> inputField;
	fieldParameters.push_back(0);
	for (auto i = 1; i <= static_cast<int>(countOfBeams); i++) {
		inputField += gaussSong2018(x1, x2, wavelength, sigma, n0, radius, velocity, relativePower, static_cast<double>(i), countOfBeams);
	}
	auto P0 = fieldParameters.at(7);
	auto k = 2 * M_PI * n0 / wavelength;
	auto gamma = sqrt(relativePower) / (sqrt(P0) * k * sigma * sigma);
	auto beta = gamma * sqrt(P0);
	auto G0 = 1 / sqrt(6 * exp((-k * k * velocity * velocity * beta * beta) / (sigma * sigma)) * (1 + 2 * exp(-3 * radius * radius * (pow(sigma, 4) - k * k * velocity * velocity * beta * beta) / (4 * sigma * sigma)) + 2 * exp(-3 * radius * radius * (pow(sigma, 4) - k * k * velocity * velocity * beta * beta) / (4 * sigma * sigma))) / pow(sigma, 4));
	return G0 * inputField;
}

vector<vector<complex<double>>> field::selectInputField(vector<double>& x1, vector<double>& x2) {
	switch (inputFunction) {
	case inputField::gauss:
		return gauss(x1, x2, fieldParameters.at(1), fieldParameters.at(2));
	case inputField::gaussHermite:
		return gaussHermite(x1, x2, fieldParameters.at(1), fieldParameters.at(2), fieldParameters.at(3));
	case inputField::gaussLaguerre:
		return gaussLaguerre(x1, x2, fieldParameters.at(1), fieldParameters.at(2), fieldParameters.at(3));
	case inputField::gaussWang2008:
		return clusterModeGaussWang2008(fieldParameters.at(0), fieldParameters.at(1), fieldParameters.at(2), fieldParameters.at(3), fieldParameters.at(4), fieldParameters.at(5));
	case inputField::gaussSong2018:
		return clusterModeGaussSong2018(x1, x2, fieldParameters.at(0), fieldParameters.at(1), fieldParameters.at(2), fieldParameters.at(3), fieldParameters.at(4), fieldParameters.at(5), fieldParameters.at(6));
	default:
		return calculatedField;
	}
}

//calculateCollinsCUDA() requires CUDA toolkit
vector<vector<complex<double>>> field::collins(vector<vector<complex<double>>>& inputFunction, vector<double>& x1, vector<double>& x2, vector<double>& x3, vector<double>& x4) {
	cout << "Используется GPU для вычисления светового поля. Пожалуйста, подождите…" << endl;
	return calculateCollinsCUDA(inputFunction, x1, x2, x3, x4, n1, n2, 2 * M_PI / fieldParameters.at(0), limits, matrixABCD);
}

vector<vector<complex<double>>> field::collinsSingular(vector<vector<complex<double>>>& inputFunction, vector<double>& x1, vector<double>& x2, vector<double>& x3, vector<double>& x4) {
	auto k = 2 * M_PI / fieldParameters.at(0);
	vector<vector<complex<double>>> output;
	for (auto p = 0; p < x4.size(); p++) {
		output.push_back(vector<complex<double>>());
		for (auto q = 0; q < x3.size(); q++) {
			output.back().push_back(matrixABCD.at(1).at(1) * inputFunction.at(p).at(q) * exp(complex<double>(0, (k * matrixABCD.at(1).at(0) * matrixABCD.at(1).at(1) * (x3.at(q) * x3.at(q) + x4.at(p) * x4.at(p))) / 2)));
		}
	}
	cout << "Сингулярный случай — световое поле посчитано мгновенно.";
	return output;
}

vector<double> field::calcPoints(double interval, double count, double shift = 0) {
	auto pointValue = -interval;
	vector<double> points;
	auto h = 2 * interval / count;
	for (auto i = 0; i < count; i++) {
		points.push_back((std::abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? pointValue * matrixABCD.at(1).at(1) - shift : pointValue - shift);
		pointValue += h;
	}
	return points;
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

double field::power(field& field) {
	auto amplitudes = abs(field);
	auto power = 0.;
	for (auto row : amplitudes) {
		for (auto value : row) {
			power += value;
		}
	}
	return power;
}

field::field() : n1(0), n2(0), limits(vector<double>()), inputFunction(inputField::undefined), selectedCrossSection(crossSection::Oxy), fieldParameters(vector<double>()), matrixABCD(vector<vector<double>>()), calculatedField(vector<vector<complex<double>>>()), superposition(false), calculated(false) {
}

field::field(vector<vector<complex<double>>>& superpositionField, crossSection crossSection) : n1(0), n2(0), limits(vector<double>()), inputFunction(inputField::undefined), selectedCrossSection(crossSection), fieldParameters(vector<double>()), matrixABCD(vector<vector<double>>()), calculatedField(superpositionField), superposition(true), calculated(true) {
}

field::field(vector<double> limits, int n1, int n2, crossSection crossSection, vector<vector<double>> matrixABCD, inputField inputFunction, vector<double> fieldParameters) : n1(n1), n2(n2), limits(limits), inputFunction(inputFunction), selectedCrossSection(crossSection), fieldParameters(fieldParameters), matrixABCD(matrixABCD), calculatedField(vector<vector<complex<double>>>()), superposition(false), calculated(false) {
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

vector<vector<complex<double>>> field::getCalculatedField() {
	return calculatedField;
}

crossSection field::getCrossSection() {
	return selectedCrossSection;
}

BMP field::createBMP(string schemeName, bool phase) {
	if (!this->isCalculated()) {
		this->calculate();
	}
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

bool field::isSuperposition() {
	return superposition;
}

bool field::isCalculated() {
	return calculated;
}

void field::calculate() {
	auto x = calcPoints(limits.at(0), n1);
	auto y = calcPoints(limits.at(1), n1);
	reverse(y.begin(), y.end());
	auto u = calcPoints(limits.at(2), n2);
	auto v = calcPoints(limits.at(3), n2);
	reverse(v.begin(), v.end());
	auto inputField = selectInputField(x, y);
	calculatedField = (std::abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? collinsSingular(inputField, x, y, u, v) : collins(inputField, x, y, u, v);
	calculated = true;
}

field operator+(field& left, field& right) {
	if (!left.isCalculated()) {
		left.calculate();
	}
	if (!right.isCalculated()) {
		right.calculate();
	}
	if ((left.getCalculatedField().at(1).size() != right.getCalculatedField().at(1).size()) || (left.getCalculatedField().size() != right.getCalculatedField().size()) || (left.getCrossSection() != right.getCrossSection())) {
		return field();
	}
	vector<vector<complex<double>>> summa;
	for (auto i = 0; i < left.getCalculatedField().at(1).size(); i++) {
		summa.push_back(vector<complex<double>>());
		for (auto j = 0; j < left.getCalculatedField().size(); j++) {
			summa.at(i).push_back(left.getCalculatedField().at(i).at(j) + right.getCalculatedField().at(i).at(j));
		}
	}
	return field(summa, left.getCrossSection());
}

vector<vector<complex<double>>> operator+(vector<vector<complex<double>>>& left, vector<vector<complex<double>>>& right) {
	if (left.size() == 0) {
		return right;
	}
	if (right.size() == 0) {
		return left;
	}
	if ((left.at(1).size() != right.at(1).size()) || (left.size() != right.size())) {
		return vector<vector<complex<double>>>();
	}
	vector<vector<complex<double>>> summa;
	for (auto i = 0; i < left.at(1).size(); i++) {
		summa.push_back(vector<complex<double>>());
		for (auto j = 0; j < left.size(); j++) {
			summa.at(i).push_back(left.at(i).at(j) + right.at(i).at(j));
		}
	}
	return summa;
}

vector<vector<complex<double>>> operator+=(vector<vector<complex<double>>>& left, vector<vector<complex<double>>>& right) {
	left = left + right;
	return left;
}

vector<vector<complex<double>>> operator*(double left, vector<vector<complex<double>>>& right) {
	for (auto i = 0; i < right.at(1).size(); i++) {
		for (auto j = 0; j < right.size(); j++) {
			right.at(i).at(j) *= left;
		}
	}
	return right;
}
