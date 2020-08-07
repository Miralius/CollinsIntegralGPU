#include "field.h"
extern vector<vector<complex<double>>> calculateCollinsCUDA(vector<vector<complex<double>>>& inputFunction, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v, int n1, int n2, double waveNumber, vector<double> limits, vector<vector<double>> matrixABCD);

vector<vector<complex<double>>> field::selectPatternField(vector<double>& x, vector<double>& y) {
	switch (pattern) {
	case patternField::solitone:
		return patternSolitone(x, y);
	case patternField::radiallySymmetricCluster:
		return patternRadiallySymmetricCluster(x, y);
	default:
		return calculatedField;
	}
}

vector<vector<complex<double>>> field::patternSolitone(vector<double>& x, vector<double>& y) {
	vector<vector<complex<double>>> input;
	for (auto i = 0; i < y.size(); i++) {
		input.push_back(vector<complex<double>>());
		for (auto j = 0; j < x.size(); j++) {
			input.back().push_back(selectInputMode(x.at(j), y.at(i), 1, i));
		}
	}
	return input;
}

vector<vector<complex<double>>> field::patternRadiallySymmetricCluster(vector<double>& x, vector<double>& y) {
	vector<vector<complex<double>>> inputField;
	auto countOfBeams = fieldParameters.at(0).at(1);
	auto radius = fieldParameters.at(0).at(2);
	auto Is_the_angle_between_the_last_beam_and_the_origin_zero = static_cast<boolean>(fieldParameters.at(0).at(3));
	for (auto k = 1; k <= static_cast<int>(countOfBeams); k++) {
		vector<vector<complex<double>>> currentMode;
		auto phi_n = Is_the_angle_between_the_last_beam_and_the_origin_zero ? M_PI * 2 * static_cast<long long>(k) / countOfBeams : M_PI * (2 * static_cast<long long>(k) - 1) / countOfBeams;
		auto c_xn = radius * cos(phi_n);
		auto c_yn = radius * sin(phi_n);
		for (auto i = 0; i < y.size(); i++) {
			currentMode.push_back(vector<complex<double>>());
			for (auto j = 0; j < x.size(); j++) {
				currentMode.back().push_back(selectInputMode(x.at(j) - c_xn, y.at(i) - c_yn, k, i));
			}
		}
		calculatedField = currentMode;
		auto P0 = power(*this);
		auto wavelength = fieldParameters.at(0).at(0);
		auto waveNumber = 2 * M_PI / wavelength;
		auto eta = fieldParameters.at(k).at(7);
		auto beamWaist = fieldParameters.at(k).at(1);
		auto gamma = sqrt(eta) / (sqrt(P0) * waveNumber * beamWaist * beamWaist);
		auto beta = gamma * sqrt(P0);
		auto ksi = fieldParameters.at(k).at(6);
		auto G0 = beamWaist * beamWaist / sqrt(6 * exp((-waveNumber * waveNumber * ksi * ksi * beta * beta) / (beamWaist * beamWaist)) * (1 + 2 * exp(-3 * radius * radius * (pow(beamWaist, 4) - waveNumber * waveNumber * ksi * ksi * beta * beta) / (4 * beamWaist * beamWaist)) + 2 * exp(-radius * radius * (pow(beamWaist, 4) - waveNumber * waveNumber * ksi * ksi * beta * beta) / (4 * beamWaist * beamWaist))));
		for (auto i = 0; i < y.size(); i++) {
			for (auto j = 0; j < x.size(); j++) {
				currentMode.at(i).at(j) *= G0 * sqrt(P0) / (sqrt(M_PI) * beamWaist) * exp(complex<double>(0, -waveNumber * ksi * beta * radius * (sin(phi_n) * x.at(j) - cos(phi_n) * y.at(i))));
			}
		}
		inputField += currentMode;
		currentMode.clear();
	}
	return inputField;
}

complex<double> field::selectInputMode(double x, double y, int beamNumber, int iterator) {
	inputField inputFunction = static_cast<inputField>(fieldParameters.at(beamNumber).at(0));
	switch (inputFunction) {
	case inputField::gauss:
		return gaussMode(x, y, beamNumber, iterator, fieldParameters.at(beamNumber).at(2), fieldParameters.at(beamNumber).at(3));
	case inputField::gaussHermite:
		return gaussHermite(x, y, beamNumber);
	case inputField::gaussLaguerre:
		return gaussLaguerre(x, y, beamNumber, iterator);
	default:
		return complex<double>(0, 0);
	}
}

complex<double> field::gaussMode(double x, double y, int beamNumber = 1, int iterator = 0, double beamWaistCoefficient = 2, double initialPhaseVariant = 0) {
	auto phi = (iterator < (n1 / 2)) ? atan2(y, -x) : (atan2(y, -x) + 2 * M_PI);
	auto countOfBeams = fieldParameters.at(0).at(1);
	auto beamWaist = fieldParameters.at(beamNumber).at(1);
	auto m = fieldParameters.at(beamNumber).at(4);
	return exp(complex<double>(-(x * x + y * y) / (beamWaistCoefficient * beamWaist * beamWaist), 0)) * (initialPhaseVariant ? exp(complex<double>(0, m * M_PI * (2 * static_cast<long long>(beamNumber) - 1) / countOfBeams)) : exp(complex<double>(0, m * phi)));
}

complex<double> field::gaussHermite(double x, double y, int beamNumber) {
	auto beamWaist = fieldParameters.at(beamNumber).at(1);
	auto m = static_cast<int>(fieldParameters.at(beamNumber).at(4));
	auto n = static_cast<int>(fieldParameters.at(beamNumber).at(5));
	//std::hermite() requires C++17
	return (hermite(n, y / beamWaist) / sqrt(pow(2, n) * tgamma(n + 1) * sqrt(M_PI))) * (hermite(m, x / beamWaist) / sqrt(pow(2, m) * tgamma(m + 1) * sqrt(M_PI))) * gaussMode(x, y);
}

complex<double> field::gaussLaguerre(double x, double y, int beamNumber, int iterator) {
	auto beamWaist = fieldParameters.at(beamNumber).at(1);
	auto m = static_cast<int>(fieldParameters.at(beamNumber).at(4));
	auto n = static_cast<int>(fieldParameters.at(beamNumber).at(5));
	//std::assoc_laguerre() requires C++17
	return (1 / beamWaist) * sqrt(2 * tgamma(n + 1) / (M_PI * tgamma(n + std::abs(m) + 1))) * gaussMode(x, y, beamNumber, iterator, 1) * pow(((sqrt(2 * (x * x + y * y))) / beamWaist), std::abs(m)) * assoc_laguerre(n, std::abs(m), 2 * (x * x + y * y) / (beamWaist * beamWaist));
}

//calculateCollinsCUDA() requires CUDA toolkit
vector<vector<complex<double>>> field::collins(vector<vector<complex<double>>>& inputFunction, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v) {
	cout << "Используется GPU для вычисления светового поля. Пожалуйста, подождите…" << endl;
	auto wavelength = fieldParameters.at(0).at(0);
	return calculateCollinsCUDA(inputFunction, x, y, u, v, n1, n2, 2 * M_PI / wavelength, limits, matrixABCD);
}

vector<vector<complex<double>>> field::collinsSingular(vector<vector<complex<double>>>& inputFunction, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v) {
	auto k = 2 * M_PI / fieldParameters.at(0).at(0);
	vector<vector<complex<double>>> output;
	for (auto p = 0; p < v.size(); p++) {
		output.push_back(vector<complex<double>>());
		for (auto q = 0; q < u.size(); q++) {
			output.back().push_back(matrixABCD.at(1).at(1) * inputFunction.at(p).at(q) * exp(complex<double>(0, (k * matrixABCD.at(1).at(0) * matrixABCD.at(1).at(1) * (u.at(q) * u.at(q) + v.at(p) * v.at(p))) / 2)));
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

void field::setFieldParameters(bool different, vector<vector<double>> parameters) {
	if (different) {
		fieldParameters = parameters;
	}
	else {
		vector<vector<double>> params;
		auto countOfBeams = parameters.at(0).at(1);
		params.push_back(parameters.at(0));
		for (auto i = 1; i <= countOfBeams; i++) {
			params.push_back(parameters.at(1));
		}
		fieldParameters = params;
	}
}

field::field() : n1(0), n2(0), limits(vector<double>()), pattern(patternField::solitone), selectedCrossSection(crossSection::Oxy), fieldParameters(vector<vector<double>>()), matrixABCD(vector<vector<double>>()), calculatedField(vector<vector<complex<double>>>()), superposition(false), calculated(false) {
}

field::field(vector<vector<complex<double>>>& superpositionField, crossSection crossSection) : n1(0), n2(0), limits(vector<double>()), pattern(patternField::solitone), selectedCrossSection(crossSection), fieldParameters(vector<vector<double>>()), matrixABCD(vector<vector<double>>()), calculatedField(superpositionField), superposition(true), calculated(true) {
}

field::field(vector<double> limits, int n1, int n2, crossSection crossSection, vector<vector<double>> matrixABCD, patternField selectedPattern) : n1(n1), n2(n2), limits(limits), pattern(selectedPattern), selectedCrossSection(crossSection), fieldParameters(vector<vector<double>>()), matrixABCD(matrixABCD), calculatedField(vector<vector<complex<double>>>()), superposition(false), calculated(false) {
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
	auto inputField = selectPatternField(x, y);
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
