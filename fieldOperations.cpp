#include "CollinsIntegralGPU.h"

vector<vector<complex<double>>> vortex(vector<vector<double>>& func, vector<double>& x, vector<double>& y, double n) {
	vector<vector<complex<double>>> functionVortex;
	for (int i = 0; i < y.size(); i++) {
		functionVortex.push_back(vector<complex<double>>());
		for (int j = 0; j < x.size(); j++) {
			functionVortex.at(i).push_back(func.at(i).at(j) * exp(complex<double>(0, (n * ((i < (y.size() / 2)) ? atan2(-y.at(i), x.at(j)) : (atan2(-y.at(i), x.at(j)) + 2 * PI))))));
		}
	}
	return functionVortex;
}

vector<vector<double>> abs(vector<vector<complex<double>>>& field) {
	vector<vector<double>> absField;
	for (vector<complex<double>> row : field) {
		absField.push_back(vector<double>());
		for (complex<double> value : row) {
			absField.back().push_back(abs(value));
		}
	}
	return absField;
}

vector<vector<double>> arg(vector<vector<complex<double>>>& field) {
	vector<vector<double>> argField;
	for (vector<complex<double>> row : field) {
		argField.push_back(vector<double>());
		for (complex<double> value : row) {
			argField.back().push_back((value.imag() < 0) ? (arg(value) + 2 * PI) : arg(value));
		}
	}
	return argField;
}

double minimum(vector<vector<double>>& field) {
	double minValue = DBL_MAX;
	for (vector<double> row : field) {
		for (double value : row) {
			minValue = min(minValue, value);
		}
	}
	return minValue;
}

double maximum(vector<vector<double>>& field) {
	double maxValue = DBL_MIN;
	for (vector<double> row : field) {
		for (double value : row) {
			maxValue = max(maxValue, value);
		}
	}
	return maxValue;
}

class pixel {
private:
	vector<unsigned char> colors;

public:
	pixel() {
		colors = vector<unsigned char>();
	}

	pixel(unsigned char blue, unsigned char green, unsigned char red, unsigned char alpha) {
		colors = { blue, green, red, alpha };
	}

	vector<unsigned char> getPixel() {
		return colors;
	}
};

istream& operator>>(istream& input, pixel& data) {
	int blue, green, red;
	input >> red >> green >> blue;
	if (!input) {
		return input;
	}
	data = pixel(blue, green, red, 255);
	return input;
}

vector<vector<unsigned char>> applyScheme(scheme schemeName) {
	vector<vector<unsigned char>> schemes;
	vector<pixel> colors;
	switch (schemeName) {
	case scheme::black_white:
		for (int i = 0; i < 256; i++) {
			schemes.push_back(vector<unsigned char>({ (unsigned char)i, (unsigned char)i, (unsigned char)i, 255 }));
		}
		break;
	case scheme::red:
		for (int i = 0; i < 256; i++) {
			schemes.push_back(vector<unsigned char>({ 0, 0, (unsigned char)i, 255 }));
		}
		break;
	case scheme::fire:
		colors = loadingData<pixel>("fire.txt");
		for (int i = 255; i >= 0; i--) {
			schemes.push_back(colors.at(i).getPixel());
		}
		break;
	default:
		error("Выбрана неверная цветовая схема!");
	}
	return schemes;
}

vector<vector<vector<unsigned char>>> fieldToBMP(vector<vector<double>> field, scheme schemeName, bool phase) {
	double minValue;
	double maxValue;
	if (phase) {
		minValue = 0;
		maxValue = 2 * PI;
	}
	else {
		minValue = 0;
		maxValue = maximum(field);
	}
	vector<vector<unsigned char>> scheme = applyScheme(schemeName);
	vector<vector<vector<unsigned char>>> pixels;
	for (vector<double> row : field) {
		pixels.push_back(vector<vector<unsigned char>>());
		for (double value : row) {
			pixels.back().push_back(scheme.at((unsigned char)round((value - minValue) * 255 / (maxValue - minValue))));
		}
	}
	return pixels;
}