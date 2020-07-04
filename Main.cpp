﻿#include "CollinsIntegralGPU.h"

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
		

		//автоматическая версия
		double a, b, c, d;
		a = 2;
		b = 2;
		c = 2;
		d = 2;
		double p = 0;
		while (p < 4.1) {
			double f(1000);
			vector<vector<double>> matrixABCD;
			matrixABCD.push_back(vector<double>());
			matrixABCD.at(0).push_back(cos(p * PI / 2));
			matrixABCD.at(0).push_back(f * sin(p * PI / 2));
			matrixABCD.push_back(vector<double>());
			matrixABCD.at(1).push_back(-sin(p * PI / 2) / f);
			matrixABCD.at(1).push_back(cos(p * PI / 2));
			int n1(500), n2(500);

			vector<double> x = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? calcPoints(c, n2, matrixABCD.at(1).at(1)) : calcPoints(a, n1);
			vector<double> y = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? calcPoints(d, n2, matrixABCD.at(1).at(1)) : calcPoints(b, n1);
			vector<double> u = calcPoints(c, n2);
			vector<double> v = calcPoints(d, n2);

			double m(5);

			vector<vector<double>> input;
			vector<vector<complex<double>>> functionVortex;
			double sigma(0.5);
			input = functionGauss(x, y, sigma);
			functionVortex = vortex(input, x, y, m);
			double wavelength(650);
			wavelength /= 1000000;
			double hx = 2 * a / n1;
			double hy = 2 * b / n1;

			vector<vector<complex<double>>> output = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? collins(functionVortex, u, v, matrixABCD, wavelength) : collins(functionVortex, x, y, u, v, matrixABCD, wavelength, hx, hy);

			BMP absInput(fieldToBMP(abs(functionVortex), scheme::fire, false));
			BMP absOutput(fieldToBMP(abs(output), scheme::fire, false));
			BMP argInput(fieldToBMP(arg(functionVortex), scheme::black_white, true));
			BMP argOutput(fieldToBMP(arg(output), scheme::black_white, true));

			writingFile<BMP>(absOutput, "absGauss_" + to_string(p) + ".bmp");
			writingFile<BMP>(argOutput, "argGauss_" + to_string(p) + ".bmp");

			cout << endl << "Создана мода Гаусса на расстоянии " << p * f << " мм." << endl;
			p += 0.5;
			matrixABCD.clear();
		}

		p = 0;
		while (p < 4.1) {
			double f(1000);
			vector<vector<double>> matrixABCD;
			matrixABCD.push_back(vector<double>());
			matrixABCD.at(0).push_back(cos(p * PI / 2));
			matrixABCD.at(0).push_back(f * sin(p * PI / 2));
			matrixABCD.push_back(vector<double>());
			matrixABCD.at(1).push_back(-sin(p * PI / 2) / f);
			matrixABCD.at(1).push_back(cos(p * PI / 2));
			int n1(500), n2(500);

			vector<double> x = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? calcPoints(c, n2, matrixABCD.at(1).at(1)) : calcPoints(a, n1);
			vector<double> y = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? calcPoints(d, n2, matrixABCD.at(1).at(1)) : calcPoints(b, n1);
			vector<double> u = calcPoints(c, n2);
			vector<double> v = calcPoints(d, n2);

			vector<vector<complex<double>>> functionVortex;
			double sigma(0.5);
			int n(0), m(1);
			int pp(2), q(0);
			functionVortex = superposition(functionGaussHermite(x, y, sigma, n, m, 0.25), functionGaussHermite(x, y, sigma, pp, q, -0.25));
			double wavelength(650);
			wavelength /= 1000000;
			double hx = 2 * a / n1;
			double hy = 2 * b / n1;

			vector<vector<complex<double>>> output = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? collins(functionVortex, u, v, matrixABCD, wavelength) : collins(functionVortex, x, y, u, v, matrixABCD, wavelength, hx, hy);

			BMP absInput(fieldToBMP(abs(functionVortex), scheme::fire, false));
			BMP absOutput(fieldToBMP(abs(output), scheme::fire, false));
			BMP argInput(fieldToBMP(arg(functionVortex), scheme::black_white, true));
			BMP argOutput(fieldToBMP(arg(output), scheme::black_white, true));

			writingFile<BMP>(absOutput, "absHermite_" + to_string(p) + ".bmp");
			writingFile<BMP>(argOutput, "argHermite_" + to_string(p) + ".bmp");
			cout << endl << "Создана суперпозиция мод Гаусса-Эрмита на расстоянии " << p * f << " мм." << endl;
			p += 0.5;
			matrixABCD.clear();
		}

		//while (1) {
		//	cout << "Введите пределы интегрирования (a, b и c, d):" << "\na = ";
		//	double a, b, c, d;
		//	while (!(cin >> a) || (a < 0)) wrongInput();
		//	if (a == 0) break;
		//	cout << "b = ";
		//	while (!(cin >> b) || !(b > 0)) wrongInput();
		//	cout << "c = ";
		//	while (!(cin >> c) || !(c > 0)) wrongInput();
		//	cout << "d = ";
		//	while (!(cin >> d) || !(d > 0)) wrongInput();
		//	cout << "Дробное преобразование Фурье (y или n)?" << "\nОтвет: ";
		//	char ff;
		//	while (!(cin >> ff) || !(ff == 'y' || ff == 'n')) wrongInput();
		//	vector<vector<double>> matrixABCD;
		//	if (ff == 'y') {
		//		double p, f;
		//		cout << "Введите дробный индекс (p) и фокусное расстояние (f, мм):" << "\np = ";
		//		while (!(cin >> p)) wrongInput();
		//		cout << "f = ";
		//		while (!(cin >> f) || !(f > 0)) wrongInput();
		//		matrixABCD.push_back(vector<double>());
		//		matrixABCD.at(0).push_back(cos(p * PI / 2));
		//		matrixABCD.at(0).push_back(f * sin(p * PI / 2));
		//		matrixABCD.push_back(vector<double>());
		//		matrixABCD.at(1).push_back(-sin(p * PI / 2) / f);
		//		matrixABCD.at(1).push_back(cos(p * PI / 2));
		//	}
		//	else {
		//		cout << "Введите коэффициенты ABCD-матрицы: " << endl;
		//		double coefficient;
		//		for (int i = 0; i < 2; i++) {
		//			matrixABCD.push_back(vector<double>());
		//			for (int j = 0; j < 2; j++) {
		//				while (!(cin >> coefficient)) wrongInput();
		//				matrixABCD.at(i).push_back(coefficient);
		//			}
		//		}
		//		if (abs((matrixABCD.at(0).at(0) * matrixABCD.at(1).at(1) - matrixABCD.at(0).at(1) * matrixABCD.at(1).at(0)) - 1) > FLT_EPSILON) {
		//			error("Определитель ABCD-матрицы должен быть равен 1");
		//		}
		//	}
		//	
		//	cout << "Введите количество отсчётов интегрирования (n входного поля и n выходного поля):" << "\nn1 = ";
		//	int n1, n2;
		//	while (!(cin >> n1) || !(n1 > 0)) wrongInput();
		//	cout << "n2 = ";
		//	while (!(cin >> n2) || !(n2 > 0)) wrongInput();

		//	vector<double> x = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? calcPoints(c, n2, matrixABCD.at(1).at(1)) : calcPoints(a, n1);
		//	vector<double> y = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? calcPoints(d, n2, matrixABCD.at(1).at(1)) : calcPoints(b, n1);
		//	vector<double> u = calcPoints(c, n2);
		//	vector<double> v = calcPoints(d, n2);

		//	cout << "Введите топологический заряд:" << "\nm = ";
		//	double m;
		//	cin >> m;
		//	
		//	vector<vector<double>> input;
		//	vector<vector<complex<double>>> functionVortex;
		//	cout << "Входная функция:" << "\nГауссов пучок (1)\nМоды Гаусса-Лагерра (n = 0, m) (2)\nСуперпозиция мод Гаусса-Эрмита (3): ";
		//	string select;
		//	cin >> select;
		//	if (select == "1") {
		//		cout << "Введите параметр сигма:" << "\nsigma = ";
		//		double sigma;
		//		while(!(cin >> sigma) || !(sigma > 0)) wrongInput();
		//		input = functionGauss(x, y, sigma);
		//		functionVortex = vortex(input, x, y, m);
		//	}
		//	else if (select == "2") {
		//		cout << "Введите параметр сигма:" << "\nsigma = ";
		//		double sigma;
		//		while (!(cin >> sigma) || !(sigma > 0)) wrongInput();
		//		input = functionGaussLaguerre(x, y, sigma, 0, m);
		//		functionVortex = vortex(input, x, y, m);
		//	}
		//	else if (select == "3") {
		//		cout << "Введите параметр сигма:" << "\nsigma = ";
		//		double sigma;
		//		while (!(cin >> sigma) || !(sigma > 0)) wrongInput();
		//		cout << "Введите порядок (n, m) первой моды:" << "\nn m = ";
		//		int n, m;
		//		while (!(cin >> n) || !(n >= 0)) wrongInput();
		//		while (!(cin >> m) || !(m >= 0)) wrongInput();
		//		cout << "Введите порядок (p, q) второй моды:" << "\np q = ";
		//		int p, q;
		//		while (!(cin >> p) || !(p >= 0)) wrongInput();
		//		while (!(cin >> q) || !(q >= 0)) wrongInput();
		//		functionVortex = superposition(functionGaussHermite(x, y, sigma, n, m, 0.25), functionGaussHermite(x, y, sigma, p, q, -0.25));
		//	}
		//	else {
		//		error("Не выбрана входная функция!");
		//	}

		//	cout << "Введите длину волны света (нм):" << "\nwavelength = ";
		//	double wavelength;
		//	cin >> wavelength;
		//	wavelength /= 1000000;
		//	double hx = 2 * a / n1;
		//	double hy = 2 * b / n1;

		//	vector<vector<complex<double>>> output = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? collins(functionVortex, u, v, matrixABCD, wavelength) : collins(functionVortex, x, y, u, v, matrixABCD, wavelength, hx, hy);
		//	//vector<vector<complex<double>>> output = (abs(matrixABCD.at(0).at(1)) < FLT_EPSILON) ? collins(functionVortex, u, v, matrixABCD, wavelength) : cuCollins(functionVortex, x, y, u, v, matrixABCD, wavelength, hx, hy);

		//	BMP absInput(fieldToBMP(abs(functionVortex), scheme::fire, false));
		//	BMP absOutput(fieldToBMP(abs(output), scheme::fire, false));
		//	BMP argInput(fieldToBMP(arg(functionVortex), scheme::black_white, true));
		//	BMP argOutput(fieldToBMP(arg(output), scheme::black_white, true));

		//	writingFile<BMP>(absInput, "absInput.bmp");
		//	writingFile<BMP>(absOutput, "absOutput.bmp");
		//	writingFile<BMP>(argInput, "argInput.bmp");
		//	writingFile<BMP>(argOutput, "argOutput.bmp");

		//	cout << endl << "Результаты записаны! Продолжить расчёты? Для выхода ввести 0" << endl;
		//}
	}
	catch (runtime_error & e) {
		cerr << "Ошибка! " << e.what() << endl;
	}
	catch (...) {
		cerr << "Неизвестная ошибка!" << endl;
	}
	return 0;
}