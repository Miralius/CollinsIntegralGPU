﻿#include <iostream>
#include <fstream>
#include <Windows.h>
#include <vector>
#include <complex>

constexpr auto PI = 3.1415926535897932384626433832795;

using namespace std;

inline void error(const string& s)
{
	throw runtime_error(s);
}

class BMP
{
	struct BMPFILEHEADER {
		unsigned short  bfType;
		unsigned long   bfSize;
		unsigned short  bfReserved1;
		unsigned short  bfReserved2;
		unsigned long   bfOffBits;
	};
	struct BMPINFOHEADER {
		unsigned long   biSize;
		long            biWidth;
		long            biHeight;
		unsigned short  biPlanes;
		unsigned short  biBitCount;
		unsigned long   biCompression;
		unsigned long   biSizeImage;
		long            biXPelsPerMeter;
		long            biYPelsPerMeter;
		unsigned long   biClrUsed;
		unsigned long   biClrImportant;
	};
	unsigned short countRGBChannel = 4;

private:
	BMPFILEHEADER bmpFileHeader;
	BMPINFOHEADER bmpInfoHeader;
	vector<vector<unsigned char>> pixels;

public:
	BMPFILEHEADER bmpFH() const { return bmpFileHeader; }
	BMPINFOHEADER bmpIH() const { return bmpInfoHeader; }
	vector<vector<unsigned char>> pxls() const { return pixels; }

	BMP() {
		bmpFileHeader = { 0, 0, 0, 0, 0 };
		bmpInfoHeader = { sizeof(BMPINFOHEADER), 0, 0, 1, (unsigned short)(countRGBChannel * 8), BI_RGB, 0, 0, 0, 0, 0 };
	}

	BMP(vector<vector<unsigned char>> picture) {
		bmpFileHeader = { 0x04D42, (unsigned long)((picture.size() * picture.size() * countRGBChannel) + sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)), 0, 0, (unsigned long)(sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)) };
		bmpInfoHeader = { sizeof(BMPINFOHEADER), (long)picture.size(), (long)picture.size(), 1, (unsigned short)(countRGBChannel * 8), BI_RGB, 0, 0, 0, 0, 0 };
		for (int i = (int)(picture.size() - 1); i >= 0; i--) {
			pixels.push_back(vector<unsigned char>());
			for (int j = 0; j < picture.size(); j++) {
				pixels.back().push_back(picture.at(i).at(j)); //blue channel
				pixels.back().push_back(picture.at(i).at(j)); //green channel
				pixels.back().push_back(picture.at(i).at(j)); //red channel
				pixels.back().push_back(255); //reserved channel
			}
		}
	}
};

template <typename T> ostream& operator<<(ostream& output, const vector<vector<T>>& data)
{
	string pixels = "";
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data.size(); j++) {
			pixels += data.at(i).at(j);
		}
	}
	return output << pixels;
}

ostream& operator<<(ostream& output, const BMP& data) {
	return output << data.bmpFH().bfType << data.bmpFH().bfSize << data.bmpFH().bfReserved1 << data.bmpFH().bfReserved2 << data.bmpFH().bfOffBits
		<< data.bmpIH().biSize << data.bmpIH().biWidth << data.bmpIH().biHeight << data.bmpIH().biPlanes << data.bmpIH().biBitCount
		<< data.bmpIH().biCompression << data.bmpIH().biSizeImage << data.bmpIH().biXPelsPerMeter << data.bmpIH().biYPelsPerMeter
		<< data.bmpIH().biClrUsed << data.bmpIH().biClrImportant << data.pxls();
}

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
	for (int i = 0; i < xy.size(); i++) {
		function1D.push_back((exp(-(xy.at(i) * xy.at(i)) / (2 * sigma * sigma))));
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
					value += functionVortex.at(x).at(y) * exp(complex<double>(0, ((k / (2 * matrixABCD.at(0).at(1))) * 
						(matrixABCD.at(0).at(0) * (xy.at(x) * xy.at(x) + xy.at(y) * xy.at(y)) - 2 * (xy.at(x) * uv.at(u) + xy.at(y) * uv.at(v)) + matrixABCD.at(1).at(1) * (uv.at(u) * uv.at(u) + uv.at(v) * uv.at(v))))));
				}
			}
			output.at(u).push_back(complex<double>(0, -(k / (2 * PI * matrixABCD.at(0).at(1)))) * value * h * h);
		}
	}
	return output;
}

vector<vector<double>> abs(vector<vector<complex<double>>> field) {
	vector<vector<double>> absField;
	for (int i = 0; i < field.size(); i++) {
		absField.push_back(vector<double>());
		for (int j = 0; j < field.size(); j++) {
			absField.at(i).push_back(abs(field.at(i).at(j)));
		}
	}
	return absField;
}

vector<vector<double>> arg(vector<vector<complex<double>>> field) {
	vector<vector<double>> argField;
	for (int i = 0; i < field.size(); i++) {
		argField.push_back(vector<double>());
		for (int j = 0; j < field.size(); j++) {
			argField.at(i).push_back(arg(field.at(i).at(j)));
		}
	}
	return argField;
}

double minimum(vector<vector<double>> field) {
	double minValue = DBL_MAX;
	for (int i = 0; i < field.size(); i++) {
		for (int j = 0; j < field.size(); j++) {
			minValue = min(minValue, field.at(i).at(j));
		}
	}
	return minValue;
}

double maximum(vector<vector<double>> field) {
	double maxValue = DBL_MIN;
	for (int i = 0; i < field.size(); i++) {
		for (int j = 0; j < field.size(); j++) {
			maxValue = max(maxValue, field.at(i).at(j));
		}
	}
	return maxValue;
}

vector<vector<unsigned char>> fieldToMonochrome(vector<vector<double>> field) {
	double minValue = minimum(field);
	double maxValue = maximum(field);
	vector<vector<unsigned char>> pixels;
	for (int i = 0; i < field.size(); i++) {
		pixels.push_back(vector<unsigned char>());
		for (int j = 0; j < field.size(); j++) {
			pixels.at(i).push_back((unsigned char)round((field.at(i).at(j) - minValue) * 255 / (maxValue - minValue)));
		}
	}
	return pixels;
}

template <typename T> void writingFile(T data, string nameFile)
{
	ofstream output(nameFile, ios::binary | ios::trunc | ios::out);
	if (!output) error("Запись в файл " + nameFile + " невозможна!");
	output << data << endl;
}

void writeFileBMP(int n) {
	HDC hdc = CreateDC(TEXT("DISPLAY"), NULL, NULL, NULL);
	HDC hdcCompatible = CreateCompatibleDC(hdc);
	DWORD dwWidth(n), dwHeight(n), dwBPP(GetDeviceCaps(hdc, BITSPIXEL)), dwNumColors(0);
	byte pBits[] = {20, 40, 60, 255, 80, 100, 120, 255, 140, 160, 180, 255, 200, 220, 240, 255};
	//HBITMAP bitmap;
	BITMAPINFO bmInfo;
	BITMAPFILEHEADER bmfh;
	bmfh.bfType = 0x04D42;
	bmfh.bfSize = ((dwWidth * dwHeight * dwBPP) / 8) + sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + (dwNumColors * sizeof(RGBQUAD));
	bmfh.bfReserved1 = 0;
	bmfh.bfReserved2 = 0;
	bmfh.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + (dwNumColors * sizeof(RGBQUAD));
	bmInfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	bmInfo.bmiHeader.biWidth = dwWidth;
	bmInfo.bmiHeader.biHeight = dwHeight;
	bmInfo.bmiHeader.biPlanes = 1;
	bmInfo.bmiHeader.biBitCount = (WORD)dwBPP;
	bmInfo.bmiHeader.biCompression = BI_RGB;
	bmInfo.bmiHeader.biSizeImage = 0;
	bmInfo.bmiHeader.biXPelsPerMeter = 0;
	bmInfo.bmiHeader.biYPelsPerMeter = 0;
	bmInfo.bmiHeader.biClrUsed = dwNumColors;
	bmInfo.bmiHeader.biClrImportant = dwNumColors;
	
	//bitmap = CreateDIBSection(hdc, &bmInfo, DIB_PAL_COLORS, &pBits, NULL, 0);
	//SetDIBits(hdc, bitmap, 0, n, &field, &bmInfo, DIB_PAL_COLORS);
	//HGDIOBJ gdiobj = SelectObject(hdcCompatible, (HGDIOBJ)bitmap);
	//BitBlt(hdcCompatible, 0, 0, dwWidth, dwHeight, hdc, 0, 0, SRCCOPY);
	
	ofstream file;
	file.open("image.bmp", ios::binary | ios::trunc | ios::out);
	file.write((char*)&bmfh, sizeof(BITMAPFILEHEADER));
	file.write((char*)&bmInfo, sizeof(BITMAPINFOHEADER));
	file.write((char*)pBits, (dwWidth * dwHeight * dwBPP) / 8);

	//DeleteObject(bitmap);
	DeleteDC(hdcCompatible);
	DeleteDC(hdc);
}

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		cout << "Расчёт одномерного интеграла Коллинза…" << endl;
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

			writingFile<BMP>(absInput, "absInput.bmp");
			writingFile<BMP>(absOutput, "absOutput.bmp");
			writingFile<BMP>(argInput, "argInput.bmp");
			writingFile<BMP>(argOutput, "argOutput.bmp");

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