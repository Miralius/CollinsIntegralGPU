#include <iostream>
#include <Windows.h>
#include <vector>
#include <complex>

constexpr auto PI = 3.1415926535897932384626433832795;

using namespace std;

inline void error(const string& s)
{
	throw runtime_error(s);
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

vector<vector<double>> functionGauss(vector<double> xy) {
	cout << "Введите параметр сигма:" << "\nsigma = ";
	double sigma;
	cin >> sigma;
	
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

void integrating2D_B0(double b, double n2, vector<vector<complex<double>>> funcVortex, double wavelength, double C, double D) {

}

void integrating2D(double b, double n2, double h, vector<vector<complex<double>>> funcVortex, vector<double> xy, double wavelength, double A, double B, double D, vector<vector<complex<double>>>& outputField) {

}

void collins2D(double a, double b, vector<vector<double>> matrixABCD, int n1, int n2, double n, double wavelength) {
	

	vector<vector<complex<double>>> outputField;
	//if (abs(matrixABCD.at(0).at(1)) < DBL_EPSILON) {
	//	integrating2D_B0(b, n2, funcVortex, wavelength, C, D);
	//} else integrating2D(b, n2, h, funcVortex, xy, wavelength, A, B, D, outputField);

	cout << "Введите имя файла результатов: ";
	string nameFile;
	cin >> nameFile;
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
			
			if ((matrixABCD.at(0).at(0) * matrixABCD.at(1).at(1) - matrixABCD.at(0).at(1) * matrixABCD.at(1).at(0)) != 1) {
				error("Определитель ABCD-матрицы должен быть равен 1");
			}

			cout << "Введите количество отсчётов интегрирования (n входного поля и n выходного поля):" << "\nn1 = ";
			int n1, n2;
			cin >> n1;
			cout << "n2 = ";
			cin >> n2;

			cout << "Введите число завихрения световой волны:" << "\nn = ";
			double n;
			cin >> n;

			cout << "Введите длину волны света:" << "\nwavelength = ";
			double wavelength;
			cin >> wavelength;

			double h = 2 * a / n1;

			vector<double> xy = (abs(matrixABCD.at(0).at(1)) < DBL_EPSILON) ? calcPoints(b, n2) : calcPoints(a, n1); //
			vector<double> uv = calcPoints(b, n2);

			vector<vector<double>> input;
			cout << "Входная функция:" << "\nГауссов пучок (1): ";
			string select;
			cin >> select;
			if (select == "1") {
				input = functionGauss(xy);
			}
			else {
				error("Не выбрана входная функция!");
			}

			vector<vector<complex<double>>> funcVortex = vortex(input, xy, n);


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
