#include "field.h"
extern vector<vector<complex<double>>> calculateCollinsCUDA(vector<vector<complex<double>>>& inputFunction, vector<double>& x1, vector<double>& x2, vector<double>& x3, vector<double>& x4, int n1, int n2, double waveNumber, vector<double> limits, vector<vector<double>> matrixABCD);

vector<vector<complex<double>>> field::gauss(vector<double>& x1, vector<double>& x2, double sigma, double m) {
	vector<vector<complex<double>>> input;
	for (auto i = 0; i < x2.size(); i++) {
		input.push_back(vector<complex<double>>());
		for (auto j = 0; j < x1.size(); j++) {
			input.back().push_back(exp(complex<double>(-(x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (2 * sigma * sigma), 0)) * exp(complex<double>(0, m * ((i < (x2.size() / 2)) ? atan2(x2.at(i), -x1.at(j)) : (atan2(x2.at(i), -x1.at(j)) + 2 * M_PI)))));
		}
	}
	return input;
}

vector<vector<complex<double>>> field::gaussFlatPhase(vector<double>& x1, vector<double>& x2, double sigma, double m, double N, double k) {
	vector<vector<complex<double>>> input;
	for (auto i = 0; i < x2.size(); i++) {
		input.push_back(vector<complex<double>>());
		for (auto j = 0; j < x1.size(); j++) {
			input.back().push_back(exp(complex<double>(0, 2 * M_PI / fieldParameters.at(0) * matrixABCD.at(0).at(1))) * exp(complex<double>(-(x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (sigma * sigma), 0)) * exp(complex<double>(0, m * M_PI * (2 * k - 1) / N)));
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
			input.back().push_back((1 / sigma) *sqrt(2 * tgamma(n + 1) / (M_PI * tgamma(n + std::abs(m) + 1))) * exp(-(x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (sigma * sigma)) * pow(((sqrt(2 * (x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)))) / sigma), std::abs(m)) * assoc_laguerre(static_cast<unsigned int>(n), static_cast<unsigned int>(std::abs(m)), 2 * (x1.at(j) * x1.at(j) + x2.at(i) * x2.at(i)) / (sigma * sigma)) * exp(complex<double>(0, m * ((i < (x2.size() / 2)) ? atan2(x2.at(i), -x1.at(j)) : (atan2(x2.at(i), -x1.at(j)) + 2 * M_PI)))));
		}
	}
	return input;
}

vector<vector<complex<double>>> field::selectInputField(vector<double>& x1, vector<double>& x2) {
	switch (inputFunction) {
	case inputField::gauss:
		return gauss(x1, x2, fieldParameters.at(1), fieldParameters.at(2));
	case inputField::gaussHermite:
		return gaussHermite(x1, x2, fieldParameters.at(1), fieldParameters.at(2), fieldParameters.at(3));
	case inputField::gaussLaguerre:
		return gaussLaguerre(x1, x2, fieldParameters.at(1), fieldParameters.at(2), fieldParameters.at(3));
	case inputField::gaussFlatPhase:
		return gaussFlatPhase(x1, x2, fieldParameters.at(1), fieldParameters.at(2), fieldParameters.at(3), fieldParameters.at(4));
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
			output.at(p).push_back(matrixABCD.at(1).at(1) * inputFunction.at(p).at(q) * exp(complex<double>(0, (k * matrixABCD.at(1).at(0) * matrixABCD.at(1).at(1) * (x3.at(q) * x3.at(q) + x4.at(p) * x4.at(p))) / 2)));
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

field::field() : n1(0), n2(0), limits(vector<double>()), inputFunction(inputField::undefined), selectedCrossSection(crossSection::Oxy), fieldParameters(vector<double>()), matrixABCD(vector<vector<double>>()), calculatedField(vector<vector<complex<double>>>()), superposition(false), calculated(false), clusterMode(false) {
}

field::field(vector<vector<complex<double>>>& superpositionField, crossSection crossSection) : n1(0), n2(0), limits(vector<double>()), inputFunction(inputField::undefined), selectedCrossSection(crossSection), fieldParameters(vector<double>()), matrixABCD(vector<vector<double>>()), calculatedField(superpositionField), superposition(true), calculated(true), clusterMode(false) {
}

field::field(vector<double> limits, int n1, int n2, crossSection crossSection, vector<vector<double>> matrixABCD, inputField inputFunction, vector<double> fieldParameters) : n1(n1), n2(n2), limits(limits), inputFunction(inputFunction), selectedCrossSection(crossSection), fieldParameters(fieldParameters), matrixABCD(matrixABCD), calculatedField(vector<vector<complex<double>>>()), superposition(false), calculated(false), clusterMode(false) {
}

vector<vector<double>> field::abs(field& field) {
	if (!field.isCalculated()) {
		field.calculate();
	}
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
	if (!field.isCalculated()) {
		field.calculate();
	}
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

bool field::isClusterMode() {
	return clusterMode;
}

void field::setClusterMode(bool selectedClusterMode) {
	clusterMode = selectedClusterMode;
}

void field::calculate() {
	vector<vector<complex<double>>> inputField;
	if (isClusterMode()) {
		auto N = static_cast<int>(fieldParameters.at(fieldParameters.size() - 2));
		auto radius = fieldParameters.back();
		for (auto i = 1; i <= N; i++) {
			fieldParameters.back() = static_cast<double>(i);
			auto x = calcPoints(limits.at(0), n1, radius * cos(M_PI * (2 * static_cast<long long>(i) - 1) / N));
			auto y = calcPoints(limits.at(1), n1, radius * sin(M_PI * (2 * static_cast<long long>(i) - 1) / N));
			reverse(y.begin(), y.end());
			inputField += selectInputField(x, y);
		}
	}
	auto x = calcPoints(limits.at(0), n1);
	auto y = calcPoints(limits.at(1), n1);
	reverse(y.begin(), y.end());
	auto u = calcPoints(limits.at(2), n2);
	auto v = calcPoints(limits.at(3), n2);
	reverse(v.begin(), v.end());
	if (!isClusterMode()) {
		inputField = selectInputField(x, y);
	}
	calculatedField = (std::abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? collinsSingular(inputField, x, y, u, v) : collins(inputField, x, y, u, v);
	calculated = true;
}

field operator+(field& firstField, field& secondField) {
	if (!firstField.isCalculated()) {
		firstField.calculate();
	}
	if (!secondField.isCalculated()) {
		secondField.calculate();
	}
	if ((firstField.getCalculatedField().at(1).size() != secondField.getCalculatedField().at(1).size()) || (firstField.getCalculatedField().size() != secondField.getCalculatedField().size()) || (firstField.getCrossSection() != secondField.getCrossSection())) {
		return field();
	}
	vector<vector<complex<double>>> summa;
	for (auto i = 0; i < firstField.getCalculatedField().at(1).size(); i++) {
		summa.push_back(vector<complex<double>>());
		for (auto j = 0; j < firstField.getCalculatedField().size(); j++) {
			summa.at(i).push_back(firstField.getCalculatedField().at(i).at(j) + secondField.getCalculatedField().at(i).at(j));
		}
	}
	return field(summa, firstField.getCrossSection());
}

vector<vector<complex<double>>> operator+(vector<vector<complex<double>>>& firstObj, vector<vector<complex<double>>>& secondObj) {
	if (firstObj.size() == 0) {
		return secondObj;
	}
	if (secondObj.size() == 0) {
		return firstObj;
	}
	if ((firstObj.at(1).size() != secondObj.at(1).size()) || (firstObj.size() != secondObj.size())) {
		return vector<vector<complex<double>>>();
	}
	vector<vector<complex<double>>> summa;
	for (auto i = 0; i < firstObj.at(1).size(); i++) {
		summa.push_back(vector<complex<double>>());
		for (auto j = 0; j < firstObj.size(); j++) {
			summa.at(i).push_back(firstObj.at(i).at(j) + secondObj.at(i).at(j));
		}
	}
	return summa;
}

vector<vector<complex<double>>> operator+=(vector<vector<complex<double>>>& firstObj, vector<vector<complex<double>>>& secondObj) {
	firstObj = firstObj + secondObj;
	return firstObj;
}
