#include <iostream>
#include <fstream>
#include <Windows.h>
#include <vector>
#include <complex>
#include <iterator>

constexpr auto PI = 3.1415926535897932384626433832795;

using namespace std;

inline void error(const string& s)
{
	throw runtime_error(s);
}

class BMP
{
	int const countRGBChannel = 4;
	int const BMPFILEHEADERsize = 14;
	int const BMPINFOHEADERsize = 40;

private:
	vector<vector<unsigned char>> pixels;
	//These vectors replace the structs. Each first element — value, each second element — size of type of var.
	vector<vector<int>> bmpFileHeader;
	vector<vector<int>> bmpInfoHeader;

	vector<unsigned char> toBinary(vector<int> number) {
		vector<unsigned char> binary;
		for (int i = 0; i < number.at(1); i++) {
			binary.push_back(number.at(0) >> (8 * i));
		}
		return binary;
	}

public:
	BMP() {
		bmpFileHeader = { {0, 2}, {0, 4}, {0, 2}, {0, 2}, {0, 4} };
		bmpInfoHeader = { {BMPINFOHEADERsize, 4}, {0, 4}, {0, 4}, {1, 2}, {countRGBChannel * 8, 2}, {BI_RGB, 4}, {0, 4}, {0, 4}, {0, 4}, {0, 4}, {0, 4} };
	}

	BMP(vector<vector<unsigned char>> picture) {
		bmpFileHeader = { {0x4D42, 2}, {(int)picture.size() * (int)picture.size() * countRGBChannel + BMPFILEHEADERsize + BMPINFOHEADERsize, 4}, {0, 2}, {0, 2}, {BMPFILEHEADERsize + BMPINFOHEADERsize, 4} };
		bmpInfoHeader = { {BMPINFOHEADERsize, 4}, {(int)picture.size(), 4}, {(int)picture.size(), 4}, {1, 2}, {countRGBChannel * 8, 2}, {BI_RGB, 4}, {0, 4}, {0, 4}, {0, 4}, {0, 4}, {0, 4} };
		for (int i = (int)(picture.size() - 1); i >= 0; i--) {
			pixels.push_back(vector<unsigned char>());
			for (unsigned char value : picture.at(i)) {
				pixels.back().push_back(value); //blue channel
				pixels.back().push_back(value); //green channel
				pixels.back().push_back(value); //red channel
				pixels.back().push_back(255); //reserved channel
			}
		}
	}

	vector<unsigned char> serialize() {
		vector<unsigned char> serializedBMP;
		serializedBMP.reserve(bmpFileHeader.at(1).at(0));
		for (vector<int> data : bmpFileHeader) {
			for (unsigned char byte : toBinary(data)) {
				serializedBMP.push_back(byte);
			}
		}
		for (vector<int> data : bmpInfoHeader) {
			for (unsigned char byte : toBinary(data)) {
				serializedBMP.push_back(byte);
			}
		}
		for (vector<unsigned char> data : pixels) {
			for (unsigned char byte : data) {
				serializedBMP.push_back(byte);
			}
		}
		return serializedBMP;
	}
};

vector<double> calcPoints(double interval, double count) {
	double pointValue = -interval;
	vector<double> points;
	double h = 2 * interval / count;

	for (int i = 0; i < count; i++) {
		points.push_back(pointValue);
		pointValue += h;
	}

	return points;
}

vector<double> calcPoints(double interval, double count, double D) {
	double pointValue = -interval;
	vector<double> points;
	double h = 2 * interval / count;

	for (int i = 0; i < count; i++) {
		points.push_back(pointValue * D);
		pointValue += h;
	}

	return points;
}

vector<vector<double>> functionGauss(vector<double> xy, double sigma) {
	vector<double> function1D;
	for (double value : xy) {
		function1D.push_back((exp(-(value * value) / (2 * sigma * sigma))));
	}

	vector<vector<double>> input;
	for (int i = 0; i < function1D.size(); i++) {
		input.push_back(vector<double>());
		for (int j = 0; j < function1D.size(); j++) {
			input.at(i).push_back(function1D.at(i) * function1D.at(j));
		}
	}

	return input;
}

vector<vector<complex<double>>> vortex(vector<vector<double>> func, vector<double> xy, double n) {
	vector<vector<complex<double>>> functionVortex;
	for (int i = 0; i < xy.size(); i++) {
		functionVortex.push_back(vector<complex<double>>());
		for (int j = 0; j < xy.size(); j++) {
			functionVortex.at(i).push_back(func.at(i).at(j) * exp(complex<double>(0, (n * ((j < (xy.size() / 2)) ? atan2(-xy.at(i), xy.at(j)) : (atan2(-xy.at(i), xy.at(j)) + 2 * PI))))));
		}
	}
	return functionVortex;
}

vector<vector<complex<double>>> collins(vector<vector<complex<double>>> functionVortex, vector<double> uv, vector<vector<double>> matrixABCD, double wavelength) {
	double k = 2 * PI / wavelength;

	vector<vector<complex<double>>> output;
	for (int u = 0; u < uv.size(); u++) {
		output.push_back(vector<complex<double>>());
		for (int v = 0; v < uv.size(); v++) {
			output.at(u).push_back(sqrt(matrixABCD.at(1).at(1)) * functionVortex.at(u).at(v) * exp(complex<double>(0, (k * matrixABCD.at(1).at(0) * matrixABCD.at(1).at(1) * (uv.at(u) + uv.at(v)) * (uv.at(u) + uv.at(v))) / 2)));
		}
	}
	return output;
}

vector<vector<complex<double>>> collins(vector<vector<complex<double>>> functionVortex, vector<double> xy, vector<double> uv, vector<vector<double>> matrixABCD, double wavelength, double h) {
	double k = 2 * PI / wavelength;
	vector<vector<complex<double>>> output;
	for (int u = 0; u < uv.size(); u++) {
		output.push_back(vector<complex<double>>());
		for (int v = 0; v < uv.size(); v++) {
			complex<double> value = 0;
			for (int x = 0; x < xy.size(); x++) {
				for (int y = 0; y < xy.size(); y++) {
					value += functionVortex.at(x).at(y) * exp(complex<double>(0, ((k / (2 * matrixABCD.at(0).at(1))) * (matrixABCD.at(0).at(0) * (xy.at(x) * xy.at(x) + xy.at(y) * xy.at(y)) - 2 * (xy.at(x) * uv.at(u) + xy.at(y) * uv.at(v)) + matrixABCD.at(1).at(1) * (uv.at(u) * uv.at(u) + uv.at(v) * uv.at(v))))));
				}
			}
			output.at(u).push_back(complex<double>(0, -(k / (2 * PI * matrixABCD.at(0).at(1)))) * value * h * h);
		}
	}
	return output;
}

vector<vector<double>> abs(vector<vector<complex<double>>> field) {
	vector<vector<double>> absField;
	for (vector<complex<double>> row : field) {
		absField.push_back(vector<double>());
		for (complex<double> value : row) {
			absField.back().push_back(abs(value));
		}
	}
	return absField;
}

vector<vector<double>> arg(vector<vector<complex<double>>> field) {
	vector<vector<double>> argField;
	for (vector<complex<double>> row : field) {
		argField.push_back(vector<double>());
		for (complex<double> value : row) {
			argField.back().push_back(arg(value));
		}
	}
	return argField;
}

double minimum(vector<vector<double>> field) {
	double minValue = DBL_MAX;
	for (vector<double> row : field) {
		for (double value : row) {
			minValue = min(minValue, value);
		}
	}
	return minValue;
}

double maximum(vector<vector<double>> field) {
	double maxValue = DBL_MIN;
	for (vector<double> row : field) {
		for (double value : row) {
			maxValue = max(maxValue, value);
		}
	}
	return maxValue;
}

vector<vector<unsigned char>> fieldToMonochrome(vector<vector<double>> field) {
	double minValue = minimum(field);
	double maxValue = maximum(field);
	vector<vector<unsigned char>> pixels;
	for (vector<double> row : field) {
		pixels.push_back(vector<unsigned char>());
		for (double value : row) {
			pixels.back().push_back((unsigned char)round((value - minValue) * 255 / (maxValue - minValue)));
		}
	}
	return pixels;
}

void writingFile(BMP picture, string nameFile)
{
	ofstream output(nameFile, ios::binary | ios::trunc | ios::out);
	vector<unsigned char> data = picture.serialize();
	if (!output) error("Запись в файл " + nameFile + " невозможна!");
	for (unsigned char value : data) {
		output << value;
	}
}

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		cout << "Расчёт двумерного интеграла Коллинза…" << endl;
		while (1) {
			cout << "Введите пределы интегрирования (a и b):" << "\na = ";
			double a, b;
			cin >> a;
			if (a == 0) break;
			cout << "b = ";
			cin >> b;

			cout << "Введите коэффициенты ABCD-матрицы: " << endl;
			double coefficient;
			vector<vector<double>> matrixABCD;
			for (int i = 0; i < 2; i++) {
				matrixABCD.push_back(vector<double>());
				for (int j = 0; j < 2; j++) {
					cin >> coefficient;
					matrixABCD.at(i).push_back(coefficient);
				}
			}
			if (abs((matrixABCD.at(0).at(0) * matrixABCD.at(1).at(1) - matrixABCD.at(0).at(1) * matrixABCD.at(1).at(0)) - 1) > DBL_EPSILON) {
				error("Определитель ABCD-матрицы должен быть равен 1");
			}

			cout << "Введите количество отсчётов интегрирования (n входного поля и n выходного поля):" << "\nn1 = ";
			int n1, n2;
			cin >> n1;
			cout << "n2 = ";
			cin >> n2;

			vector<double> xy = (abs(matrixABCD.at(0).at(1)) < DBL_EPSILON) ? calcPoints(b, n2, matrixABCD.at(1).at(1)) : calcPoints(a, n1);
			vector<double> uv = calcPoints(b, n2);

			vector<vector<double>> input;
			cout << "Входная функция:" << "\nГауссов пучок (1): ";
			string select;
			cin >> select;
			if (select == "1") {
				cout << "Введите параметр сигма:" << "\nsigma = ";
				double sigma;
				cin >> sigma;
				input = functionGauss(xy, sigma);
			}
			else {
				error("Не выбрана входная функция!");
			}

			cout << "Введите число завихрения световой волны:" << "\nn = ";
			double n;
			cin >> n;
			vector<vector<complex<double>>> functionVortex = vortex(input, xy, n);

			cout << "Введите длину волны света:" << "\nwavelength = ";
			double wavelength;
			cin >> wavelength;
			double h = 2 * a / n1;
			vector<vector<complex<double>>> output = (abs(matrixABCD.at(0).at(1)) < DBL_EPSILON) ? collins(functionVortex, uv, matrixABCD, wavelength) : collins(functionVortex, xy, uv, matrixABCD, wavelength, h);

			BMP absInput(fieldToMonochrome(abs(functionVortex)));
			BMP absOutput(fieldToMonochrome(abs(output)));
			BMP argInput(fieldToMonochrome(arg(functionVortex)));
			BMP argOutput(fieldToMonochrome(arg(output)));

			writingFile(absInput, "absInput.bmp");
			writingFile(absOutput, "absOutput.bmp");
			writingFile(argInput, "argInput.bmp");
			writingFile(argOutput, "argOutput.bmp");

			cout << "Продолжить расчёты? Для выхода ввести 0" << endl;
		}
	}
	catch (runtime_error & e) {
		cerr << "Ошибка! " << e.what() << endl;
	}
	catch (...) {
		cerr << "Неизвестная ошибка!" << endl;
	}
	return 0;
}