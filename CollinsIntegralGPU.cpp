#include "CollinsIntegralGPU.h"

void processing(int NOW, int MAX, int seconds, int timeLeft) {
	float proc, nowf, maxf;
	if (NOW == MAX) proc = 100.;
	else {
		nowf = (float)NOW;
		maxf = (float)MAX;
		proc = trunc(10000 * (nowf / maxf)) / 100;
	}
	cout << '\r' << "Выполнено " << setw(6) << proc << "%, прошло " << setw(6) << seconds << " секунд, осталось " << setw(6) << timeLeft << " секунд";
}

class BMP {
	int const countRGBChannel = 4;
	int const BMPFILEHEADERsize = 14;
	int const BMPINFOHEADERsize = 40;

private:
	vector<vector<vector<unsigned char>>> pixels;
	//These vectors replace the structs. Each first element — value, each second element — size of type of var.
	vector<vector<int>> bmpFileHeader;
	vector<vector<int>> bmpInfoHeader;

	void initHeaders(int width, int height) {
		bmpFileHeader = { {0x4D42, 2}, {width * height * countRGBChannel + BMPFILEHEADERsize + BMPINFOHEADERsize, 4}, {0, 2}, {0, 2}, {BMPFILEHEADERsize + BMPINFOHEADERsize, 4} };
		bmpInfoHeader = { {BMPINFOHEADERsize, 4}, {width, 4}, {height, 4}, {1, 2}, {countRGBChannel * 8, 2}, {BI_RGB, 4}, {0, 4}, {0, 4}, {0, 4}, {0, 4}, {0, 4} };
	}

	vector<unsigned char> toBinary(vector<int> number) {
		vector<unsigned char> binary;
		for (int i = 0; i < number.at(1); i++) {
			binary.push_back(number.at(0) >> (8 * i));
		}
		return binary;
	}

	vector<int> toNumber(vector<unsigned char> binary) {
		vector<int> number;
		int temp = 0;
		for (int i = 0; i < binary.size(); i++) {
			temp |= binary.at(i) << (8 * i);
		}
		number.push_back(temp);
		number.push_back(binary.size());
		return number;
	}

public:
	BMP() {
		initHeaders(0, 0);
	}

	BMP(vector<vector<vector<unsigned char>>> picture) : pixels(picture) {
		initHeaders((int)picture.size(), (int)picture.at(0).size());
	}

	//BMP(string filename, int width, int height) : BMP() { //repair this too
	//	unsigned char buf;
	//	vector<unsigned char> buf_vector;
	//	ifstream input(filename, ios::binary | ios::in);
	//	if (!input) error("Считывание с " + filename + " невозможно!");
	//	for (vector<int> data : bmpFileHeader) {
	//		for (int i = 0; i < data.at(1); i++) {
	//			input >> buf;
	//			buf_vector.push_back(buf);
	//		}
	//		data = toNumber(buf_vector);
	//		buf_vector.clear();
	//	}
	//	buf_vector.clear();
	//	for (vector<int> data : bmpInfoHeader) {
	//		for (int i = 0; i < data.at(1); i++) {
	//			input >> buf;
	//			buf_vector.push_back(buf);
	//		}
	//		data = toNumber(buf_vector);
	//		buf_vector.clear();
	//	}
	//	buf_vector.clear();
	//	for (int i = 0; i < height; i++) {
	//		pixels.push_back(vector<unsigned char>());
	//		for (int j = 0; j < width; j++) {
	//			input >> buf;
	//			pixels.at(i).push_back(255); //reserved channel
	//			input >> buf;
	//			pixels.at(i).push_back(buf); //red channel
	//			input >> buf;
	//			pixels.at(i).push_back(buf); //green channel
	//			input >> buf;
	//			pixels.at(i).push_back(buf); //blue channel
	//		}
	//	}
	//	reverse(pixels.begin(), pixels.end());
	//}

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
		for_each(pixels.rbegin(), pixels.rend(), [&](vector<vector<unsigned char>> row) {
			for (vector<unsigned char> pixel : row) {
				for (unsigned char color : pixel) {
					serializedBMP.push_back(color);
				}
			}
		});
		return serializedBMP;
	}

};

ostream& operator<<(ostream& output, BMP& bmp)
{
	vector<unsigned char> data = bmp.serialize();
	for (unsigned char value : data) {
		output << value;
	}
	return output;
}

enum class scheme {
	black_white, red
};

vector<vector<unsigned char>> applyScheme(scheme schemeName) {
	vector<vector<unsigned char>> schemes;
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
	default:
		error("Выбрана неверная цветовая схема!");
	}
	return schemes;
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

vector<vector<double>> functionGauss(vector<double> x, vector<double> y, double sigma) {
	vector<vector<double>> input;
	for (int i = 0; i < y.size(); i++) {
		input.push_back(vector<double>());
		for (int j = 0; j < x.size(); j++) {
			input.at(i).push_back((exp(-(x.at(j) * x.at(j) + y.at(i) * y.at(i)) / (2 * sigma * sigma))));
		}
	}
	return input;
}

vector<vector<double>> functionGaussLaguerre(vector<double> x, vector<double> y, double sigma, int n, double m) {
	vector<vector<double>> input;
	for (int i = 0; i < y.size(); i++) {
		input.push_back(vector<double>());
		for (int j = 0; j < x.size(); j++) {
			input.at(i).push_back((exp(-(x.at(j) * x.at(j) + y.at(i) * y.at(i)) / (2 * sigma * sigma))) * pow(((sqrt(x.at(j) * x.at(j) + y.at(i) * y.at(i))) / sigma), (int) abs(m)));
		}
	}
	return input;
}

vector<vector<complex<double>>> vortex(vector<vector<double>> func, vector<double> x, vector<double> y, double n) {
	vector<vector<complex<double>>> functionVortex;
	for (int i = 0; i < y.size(); i++) {
		functionVortex.push_back(vector<complex<double>>());
		for (int j = 0; j < x.size(); j++) {
			functionVortex.at(i).push_back(func.at(i).at(j) * exp(complex<double>(0, (n * ((i < (y.size() / 2)) ? atan2(-y.at(i), x.at(j)) : (atan2(-y.at(i), x.at(j)) + 2 * PI))))));
		}
	}
	return functionVortex;
}

vector<vector<complex<double>>> collins(vector<vector<complex<double>>> functionVortex, vector<double> u, vector<double> v, vector<vector<double>> matrixABCD, double wavelength) {
	double k = 2 * PI / wavelength;

	vector<vector<complex<double>>> output;
	for (int i = 0; i < v.size(); i++) {
		output.push_back(vector<complex<double>>());
		for (int j = 0; j < u.size(); j++) {
			output.at(i).push_back(sqrt(matrixABCD.at(1).at(1)) * functionVortex.at(i).at(j) * exp(complex<double>(0, (k * matrixABCD.at(1).at(0) * matrixABCD.at(1).at(1) * (u.at(j) * u.at(j) + v.at(i) + v.at(i))) / 2)));
		}
	}
	return output;
}

vector<vector<complex<double>>> collins(vector<vector<complex<double>>> functionVortex, vector<double> x, vector<double> y, vector<double> u, vector<double> v, vector<vector<double>> matrixABCD, double wavelength, double hx, double hy) {
	double k = 2 * PI / wavelength;
	int startTime(clock() / CLOCKS_PER_SEC), endTime(clock() / CLOCKS_PER_SEC), currentTime(clock() / CLOCKS_PER_SEC), progress(0);
	vector<vector<complex<double>>> output;
	for (int p = 0; p < v.size(); p++) {
		processing(++progress, (int)v.size(), clock() / CLOCKS_PER_SEC - currentTime, (endTime - startTime) * ((int)v.size() - p));
		startTime = endTime;
		output.push_back(vector<complex<double>>());
		for (int q = 0; q < u.size(); q++) {
			complex<double> value = 0;
			for (int i = 0; i < y.size(); i++) {
				for (int j = 0; j < x.size(); j++) {
					value += functionVortex.at(i).at(j) * exp(complex<double>(0, ((k / (2 * matrixABCD.at(0).at(1))) * (matrixABCD.at(0).at(0) * (y.at(i) * y.at(i) + x.at(j) * x.at(j)) - 2 * (y.at(i) * v.at(p) + x.at(j) * u.at(q)) + matrixABCD.at(1).at(1) * (v.at(p) * v.at(p) + u.at(q) * u.at(q))))));
				}
			}
			output.at(p).push_back(complex<double>(0, -(k / (2 * PI * matrixABCD.at(0).at(1)))) * value * hx * hy);
		}
		endTime = clock() / CLOCKS_PER_SEC;
	}
	return output;
}

vector<vector<complex<double>>> cuCollins(vector<vector<complex<double>>> functionVortex, vector<double> x, vector<double> y, vector<double> u, vector<double> v, vector<vector<double>> matrixABCD, double wavelength, double hx, double hy) {
#pragma warning(push)
#pragma warning(disable:6386)
	double k = 2 * PI / wavelength;
	
	complex<double>** functionVortexTemp = new complex<double>* [functionVortex.size() * functionVortex.at(0).size()];
	for (int i = 0; i < functionVortex.size(); i++) {
		functionVortexTemp[i] = new complex<double>[functionVortex.size()];
		for (int j = 0; j < functionVortex.at(i).size(); j++) {
			functionVortexTemp[i][j] = functionVortex.at(i).at(j);
		}
	}

	double* xTemp = new double[x.size()];
	for (int i = 0; i < x.size(); i++) {
		xTemp[i] = x.at(i);
	}

	double* yTemp = new double[y.size()];
	for (int i = 0; i < y.size(); i++) {
		yTemp[i] = y.at(i);
	}

	double* uTemp = new double[u.size()];
	for (int i = 0; i < u.size(); i++) {
		uTemp[i] = u.at(i);
	}

	double* vTemp = new double[v.size()];
	for (int i = 0; i < v.size(); i++) {
		vTemp[i] = v.at(i);
	}

	double* constantStorage = new double[8];
	constantStorage[0] = matrixABCD.at(0).at(0);
	constantStorage[1] = matrixABCD.at(0).at(1);
	constantStorage[2] = matrixABCD.at(1).at(0);
	constantStorage[3] = matrixABCD.at(0).at(1);
	constantStorage[4] = k;
	constantStorage[5] = wavelength;
	constantStorage[6] = hx;
	constantStorage[7] = hy;


#pragma warning(pop)
	vector<vector<complex<double>>> output;
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
			argField.back().push_back((value.imag() < 0) ? (arg(value) + 2 * PI) : arg(value));
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

vector<vector<vector<unsigned char>>> fieldToBMP(vector<vector<double>> field, scheme schemeName) {
	double minValue = minimum(field);
	double maxValue = maximum(field);
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

void writingFile(BMP picture, string nameFile) {
	ofstream output(nameFile, ios::binary | ios::trunc | ios::out);
	//vector<unsigned char> data = picture.serialize();
	if (!output) error("Запись в файл " + nameFile + " невозможна!");
	//for (unsigned char value : data) {
	//	output << value;
	//}
	output << picture;
}

void wrongInput() {
	cout << "Неверный ввод! Введите ещё раз!" << endl;
	cin.clear();
	cin.ignore(cin.rdbuf()->in_avail(), '\n');
}

//__global__
//void add(int n, float* x, float* y)
//{
//	int index = blockIdx.x * blockDim.x + threadIdx.x;
//	int stride = blockDim.x * gridDim.x;
//	for (int i = index; i < n; i += stride)
//		y[i] = x[i] + y[i];
//}

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		cout << "Расчёт двумерного интеграла Коллинза…" << endl;

		//BMP test("argInput.bmp", 100, 100);

		//int N = 1 << 20; // 1M elements

		//float* x, * y;
		//cudaMallocManaged(&x, N * sizeof(float));
		//cudaMallocManaged(&y, N * sizeof(float));

		//// initialize x and y arrays on the host
		//for (int i = 0; i < N; i++) {
		//	x[i] = 1.0f;
		//	y[i] = 2.0f;
		//}

		//// Run kernel on 1M elements on the CPU
		//add<<<1, 1>>>(N, x, y);

		//cudaDeviceSynchronize();

		//// Check for errors (all values should be 3.0f)
		//float maxError = 0.0f;
		//for (int i = 0; i < N; i++)
		//	maxError = fmax(maxError, fabs(y[i] - 3.0f));
		//std::cout << "Max error: " << maxError << std::endl;

		//// Free memory
		//cudaFree(x);
		//cudaFree(y);


		while (1) {
			cout << "Введите пределы интегрирования (a, b и c, d):" << "\na = ";
			double a, b, c, d;
			while (!(cin >> a) || (a < 0)) wrongInput();
			if (a == 0) break;
			cout << "b = ";
			while (!(cin >> b) || !(b > 0)) wrongInput();
			cout << "c = ";
			while (!(cin >> c) || !(c > 0)) wrongInput();
			cout << "d = ";
			while (!(cin >> d) || !(d > 0)) wrongInput();
			cout << "Дробное преобразование Фурье (y или n)?" << "\nОтвет: ";
			char ff;
			while (!(cin >> ff) || !(ff == 'y' || ff == 'n')) wrongInput();
			vector<vector<double>> matrixABCD;
			if (ff == 'y') {
				double p, f;
				cout << "Введите дробный индекс (p) и фокусное расстояние (f, мм):" << "\np = ";
				while (!(cin >> p)) wrongInput();
				cout << "f = ";
				while (!(cin >> f) || !(f > 0)) wrongInput();
				matrixABCD.push_back(vector<double>());
				matrixABCD.at(0).push_back(cos(p * PI / 2));
				matrixABCD.at(0).push_back(f * sin(p * PI / 2));
				matrixABCD.push_back(vector<double>());
				matrixABCD.at(1).push_back(-sin(p * PI / 2) / f);
				matrixABCD.at(1).push_back(cos(p * PI / 2));
			}
			else {
				cout << "Введите коэффициенты ABCD-матрицы: " << endl;
				double coefficient;
				for (int i = 0; i < 2; i++) {
					matrixABCD.push_back(vector<double>());
					for (int j = 0; j < 2; j++) {
						while (!(cin >> coefficient)) wrongInput();
						matrixABCD.at(i).push_back(coefficient);
					}
				}
				if (abs((matrixABCD.at(0).at(0) * matrixABCD.at(1).at(1) - matrixABCD.at(0).at(1) * matrixABCD.at(1).at(0)) - 1) > DBL_EPSILON) {
					error("Определитель ABCD-матрицы должен быть равен 1");
				}
			}
			
			cout << "Введите количество отсчётов интегрирования (n входного поля и n выходного поля):" << "\nn1 = ";
			int n1, n2;
			while (!(cin >> n1) || !(n1 > 0)) wrongInput();
			cout << "n2 = ";
			while (!(cin >> n2) || !(n2 > 0)) wrongInput();

			vector<double> x = (abs(matrixABCD.at(0).at(1)) < DBL_EPSILON) ? calcPoints(c, n2, matrixABCD.at(1).at(1)) : calcPoints(a, n1);
			vector<double> y = (abs(matrixABCD.at(0).at(1)) < DBL_EPSILON) ? calcPoints(d, n2, matrixABCD.at(1).at(1)) : calcPoints(b, n1);
			vector<double> u = calcPoints(c, n2);
			vector<double> v = calcPoints(d, n2);

			cout << "Введите топологический заряд:" << "\nm = ";
			double m;
			cin >> m;
			
			vector<vector<double>> input;
			cout << "Входная функция:" << "\nГауссов пучок (1)\nМоды Гаусса-Лагерра (n = 0, m) (2): ";
			string select;
			cin >> select;
			if (select == "1") {
				cout << "Введите параметр сигма:" << "\nsigma = ";
				double sigma;
				while(!(cin >> sigma) || !(sigma > 0)) wrongInput();
				input = functionGauss(x, y, sigma);
			}
			else if (select == "2") {
				cout << "Введите параметр сигма:" << "\nsigma = ";
				double sigma;
				while (!(cin >> sigma) || !(sigma > 0)) wrongInput();
				input = functionGaussLaguerre(x, y, sigma, 0, m);
			}
			else {
				error("Не выбрана входная функция!");
			}

			vector<vector<complex<double>>> functionVortex = vortex(input, x, y, m);

			cout << "Введите длину волны света (нм):" << "\nwavelength = ";
			double wavelength;
			cin >> wavelength;
			wavelength /= 1000000;
			double hx = 2 * a / n1;
			double hy = 2 * b / n1;

			vector<vector<complex<double>>> output = (abs(matrixABCD.at(0).at(1)) < DBL_EPSILON) ? collins(functionVortex, u, v, matrixABCD, wavelength) : collins(functionVortex, x, y, u, v, matrixABCD, wavelength, hx, hy);
			//vector<vector<complex<double>>> output = (abs(matrixABCD.at(0).at(1)) < DBL_EPSILON) ? collins(functionVortex, u, v, matrixABCD, wavelength) : cuCollins(functionVortex, x, y, u, v, matrixABCD, wavelength, hx, hy);

			BMP absInput(fieldToBMP(abs(functionVortex), scheme::red));
			BMP absOutput(fieldToBMP(abs(output), scheme::red));
			BMP argInput(fieldToBMP(arg(functionVortex), scheme::black_white));
			BMP argOutput(fieldToBMP(arg(output), scheme::black_white));

			writingFile(absInput, "absInput.bmp");
			writingFile(absOutput, "absOutput.bmp");
			writingFile(argInput, "argInput.bmp");
			writingFile(argOutput, "argOutput.bmp");

			cout << endl << "Результаты записаны! Продолжить расчёты? Для выхода ввести 0" << endl;
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