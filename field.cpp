#include "field.h"
extern vector<vector<complex<double>>> calculateCollinsCUDA(vector<vector<complex<double>>>& inputFunction, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v, int n1, int n2, double waveNumber, vector<double> limits, vector<vector<double>> matrixABCD);

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

vector<vector<complex<double>>> field::shift(double x0, double y0) {
	vector<vector<complex<double>>> shiftedField;
	auto oldX = calcPoints(limits.at(0), n1);
	auto oldY = calcPoints(limits.at(1), n1);
	reverse(oldY.begin(), oldY.end());
	auto newX = calcPoints(limits.at(0), n1, x0);
	auto newY = calcPoints(limits.at(1), n1, y0);
	reverse(newY.begin(), newY.end());
	int leftX;
	int rightX;
	int downY;
	int upY;
	if (x0 > 0) {
		leftX = 0;
		rightX = 0;
		do {
			rightX++;
		} while (oldX.at(rightX) < (limits.at(0) - x0));
	}
	else {
		rightX = n1;
		leftX = 0;
		while (oldX.at(leftX) < (-limits.at(0) + x0)) {
			leftX++;
		}
	}
	if (y0 > 0) {
		downY = n1;
		upY = 0;
		while (oldY.at(upY) > (limits.at(1) - y0)) {
			upY++;
		}
	}
	else {
		upY = 0;
		downY = 0;
		do {
			downY++;
		} while (oldY.at(downY) > (-limits.at(1) + y0));
	}
	for (auto i = 0; i < n1; i++) {
		shiftedField.push_back(vector<complex<double>>());
		for (auto j = 0; j < n1; j++) {
			if (x0 > 0) {
				if (y0 > 0) {
					shiftedField.back().push_back((((i + (n1 - rightX)) <= n1) & ((j + upY) <= n1)) ? calculatedField.at(i + (n1 - rightX)).at(j + upY) : 0);
				}
				else {
					shiftedField.back().push_back((((i + (n1 - rightX)) <= n1) & ((j + upY) <= n1)) ? calculatedField.at(i + (n1 - rightX)).at(j + upY) : 0);
				}
			}
			
		}
	}
	return shiftedField;
}

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
			input.back().push_back(getNormalizedCoefficient(1) * selectInputMode(x.at(j), y.at(i), 1, i));
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
				currentMode.back().push_back((static_cast<int>(fieldParameters.at(k).at(0)) == 6) ? selectInputMode(x.at(j), y.at(i), k, i) : selectInputMode(x.at(j) - c_xn, y.at(i) - c_yn, k, i));
			}
		}
		calculatedField = currentMode;
		if (static_cast<int>(fieldParameters.at(k).at(0)) == 6) {
			auto matrixABCDtemp = matrixABCD;
			matrixABCD = { {0., 1000.}, {-0.001, 0} };
			auto inputFieldTemp = calculatedField;
			calculatedField = (std::abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? collinsSingular(inputFieldTemp, x, y, x, y) : collins(inputFieldTemp, x, y, x, y);
			matrixABCD = matrixABCDtemp;
			currentMode = calculatedField;
		}
		auto P0 = power(*this);
		auto wavelength = fieldParameters.at(0).at(0);
		auto waveNumber = 2 * M_PI / wavelength;
		auto eta = fieldParameters.at(k).at(9);
		auto beamWaist = fieldParameters.at(k).at(1);
		auto gamma = sqrt(eta) / (sqrt(P0) * waveNumber * beamWaist * beamWaist);
		auto beta = gamma * sqrt(P0);
		auto ksi = fieldParameters.at(k).at(8);
		auto normalizedCoefficient = getNormalizedCoefficient(k, P0);
		for (auto i = 0; i < y.size(); i++) {
			for (auto j = 0; j < x.size(); j++) {
				currentMode.at(i).at(j) *= normalizedCoefficient * sqrt(P0) / (sqrt(M_PI) * beamWaist) * exp(complex<double>(0, -waveNumber * ksi * beta * radius * (sin(phi_n) * x.at(j) - cos(phi_n) * y.at(i))));
				currentMode.at(i).at(j) *= (static_cast<int>(fieldParameters.at(k).at(0)) == 6) ? gaussMode(x.at(j), y.at(i), k, 0) : 1;
			}
		}
		inputField += currentMode;
		currentMode.clear();
	}
	return inputField;
}

complex<double> field::getNormalizedCoefficient(int beamNumber, double inputPower) {
	auto normalizedCoefficientVariant = static_cast<int>(fieldParameters.at(beamNumber).at(10));
	switch (normalizedCoefficientVariant) {
	case 1:
		return 1;
	case 2:
		return getNormalizedCoefficientWang2008();
	case 3:
		return getNormalizedCoefficientSong2018(beamNumber, inputPower);
	case 4:
		return getNormalizedCoefficientSong2019(beamNumber, inputPower);
	default:
		return complex<double>();
	}
}

complex<double> field::getNormalizedCoefficientWang2008() {
	auto k = 2 * M_PI / fieldParameters.at(0).at(0);
	auto L0 = fieldParameters.at(0).at(4);
	return exp(complex<double>(0, k * L0));
}

complex<double> field::getNormalizedCoefficientSong2018(int beamNumber, double inputPower) {
	auto wavelength = fieldParameters.at(0).at(0);
	auto radius = fieldParameters.at(0).at(2);
	auto waveNumber = 2 * M_PI / wavelength;
	auto eta = fieldParameters.at(beamNumber).at(9);
	auto beamWaist = fieldParameters.at(beamNumber).at(1);
	auto gamma = sqrt(eta) / (sqrt(inputPower) * waveNumber * beamWaist * beamWaist);
	auto beta = gamma * sqrt(inputPower);
	auto ksi = fieldParameters.at(beamNumber).at(8);
	return beamWaist * beamWaist / sqrt(6 * exp((-waveNumber * waveNumber * ksi * ksi * beta * beta) / (beamWaist * beamWaist)) * (1 + 2 * exp(-3 * radius * radius * (pow(beamWaist, 4) - waveNumber * waveNumber * ksi * ksi * beta * beta) / (4 * beamWaist * beamWaist)) + 2 * exp(-radius * radius * (pow(beamWaist, 4) - waveNumber * waveNumber * ksi * ksi * beta * beta) / (4 * beamWaist * beamWaist))));
}

complex<double> field::getNormalizedCoefficientSong2019(int beamNumber, double inputPower) {
	auto beamWaist = fieldParameters.at(beamNumber).at(1);
	auto m = static_cast<int>(fieldParameters.at(beamNumber).at(4));
	auto n = static_cast<int>(fieldParameters.at(beamNumber).at(5));
	return sqrt(pow(2, 2 * n + m + 1) * inputPower / (M_PI * beamWaist * beamWaist * tgamma(2 * n + m + 1)));
}

complex<double> field::selectInputMode(double x, double y, int beamNumber, int iterator) {
	inputField inputFunction = static_cast<inputField>(fieldParameters.at(beamNumber).at(0));
	switch (inputFunction) {
	case inputField::gauss:
		return gaussMode(x, y, beamNumber, iterator);
	case inputField::gaussHermite:
		return gaussHermite(x, y, beamNumber);
	case inputField::gaussLaguerre:
		return gaussLaguerre(x, y, beamNumber, iterator);
	case inputField::AVB:
		return AVB(x, y, beamNumber, iterator);
	case inputField::airy:
		return airyMode(x, y, beamNumber);
	case inputField::airyShift:
		return airyShiftMode(x, y, beamNumber);
	default:
		return complex<double>(0, 0);
	}
}

complex<double> field::gaussMode(double x, double y, int beamNumber, int iterator) {
	auto phi = (iterator < (n1 / 2)) ? atan2(y, -x) : (atan2(y, -x) + 2 * M_PI);
	auto countOfBeams = fieldParameters.at(0).at(1);
	auto beamWaist = fieldParameters.at(beamNumber).at(1);
	auto beamWaistCoefficient = fieldParameters.at(beamNumber).at(2);
	auto initialPhaseVariant = static_cast<int>(fieldParameters.at(beamNumber).at(3));
	auto m = fieldParameters.at(beamNumber).at(4);
	complex<double> initialPhase;
	switch (initialPhaseVariant) {
	case 0:
		initialPhase = exp(complex<double>(0, m * phi));
		break;
	case 1:
		initialPhase = exp(complex<double>(0, m * M_PI * (2 * static_cast<long long>(beamNumber) - 1) / countOfBeams));
		break;
	case 2:
		initialPhase = exp(complex<double>(0, -m * phi));
		break;
	default:
		initialPhase = 0;
	}
	return exp(complex<double>(-(x * x + y * y) / (beamWaistCoefficient * beamWaist * beamWaist), 0)) * initialPhase;
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
	return (1 / beamWaist) * sqrt(2 * tgamma(n + 1) / (M_PI * tgamma(n + std::abs(m) + 1))) * gaussMode(x, y, beamNumber, iterator) * pow(((sqrt(2 * (x * x + y * y))) / beamWaist), std::abs(m)) * assoc_laguerre(n, std::abs(m), 2 * (x * x + y * y) / (beamWaist * beamWaist));
	//return gaussMode(x, y, beamNumber, iterator) * pow(((sqrt((x * x + y * y)))), std::abs(m)) * assoc_laguerre(n, std::abs(m), (x * x + y * y));
}

complex<double> field::AVB(double x, double y, int beamNumber, int iterator) {
	auto beamWaist = fieldParameters.at(beamNumber).at(1);
	auto m = static_cast<int>(fieldParameters.at(beamNumber).at(4));
	auto n = static_cast<int>(fieldParameters.at(beamNumber).at(5));
	return pow(sqrt(x * x + y * y) / beamWaist, 2 * n + std::abs(m)) * gaussMode(x, y, beamNumber, iterator);
}

complex<double> field::airyMode(double x, double y, int beamNumber) {
	/*auto airy = [](double x) {
		if (x > 0) {
			return M_1_PI * sqrt(x / 3) * cyl_bessel_k(1. / 3, (2. / 3) * pow(x, 3. / 2));
		}
		else if (x < 0) {
			return (sqrt(-x) / 3) * (cyl_bessel_j(1. / 3, (2. / 3) * pow(-x, 3. / 2)) + cyl_bessel_j(-1. / 3, (2. / 3) * pow(-x, 3. / 2)));
		}
		else {
			return 1 / (pow(3, 2. / 3) * tgamma(2. / 3));
		}
	};
	return airy(x) * airy(y) * gaussMode(x, y, beamNumber, 0);*/
	auto a = [](double phi) {
		return exp(complex<double>(0, phi));
	};
	complex<double> integ = 0;
	auto h = 2 * M_PI / n1;
	double phi = 0;
	while (phi <= 2 * M_PI) {
		integ += a(phi) * exp(complex<double>(0, beamNumber * (x * cos(phi) + y * sin(phi))));
		phi += h;
	}
	auto k = 2 * M_PI / fieldParameters.at(0).at(0);
	auto k_z = sqrt(k * k - beamNumber * beamNumber);
	return integ * exp(complex<double>(0, k_z * 1000));
}

complex<double> field::airyShiftMode(double x, double y, int beamNumber) {
	auto alpha = fieldParameters.at(beamNumber).at(4);
	auto beta = fieldParameters.at(beamNumber).at(5);
	auto alpha0 = fieldParameters.at(beamNumber).at(6);
	auto beta0 = fieldParameters.at(beamNumber).at(7);
	return exp(complex<double>(0, alpha * pow(x, 3) + beta * pow(y, 3) + alpha0 * x + beta0 * y));
}

//calculateCollinsCUDA() requires CUDA toolkit
vector<vector<complex<double>>> field::collins(vector<vector<complex<double>>>& inputFunction, vector<double>& x, vector<double>& y, vector<double>& u, vector<double>& v) {
	cout << "Используется GPU для вычисления светового поля. Пожалуйста, подождите…" << endl;
	auto wavelength = fieldParameters.at(0).at(0);
	auto field = calculateCollinsCUDA(inputFunction, x, y, u, v, n1, n2, 2 * M_PI / wavelength, limits, matrixABCD);
	cout << "\rВыполнено 100.00%";
	return field;
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
