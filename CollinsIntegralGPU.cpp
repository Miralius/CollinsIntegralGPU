#include <iostream>
#include <Windows.h>
#include <vector>

constexpr auto PI = 3.1415926535897932384626433832795;

using namespace std;

inline void error(const string& s)
{
	throw runtime_error(s);
}

class Complex {
private:
	double value[2];
	
public:
	double real() const {
		return value[0];
	}
	
	double imag() const { 
		return value[1];
	}

	Complex() {
		value[0] = 0;
		value[1] = 0;
	}
	
	Complex& operator=(const double right) {
		value[0] = right;
		value[1] = 0;
		return *this;
	}
};

void vortex(vector<vector<double>> func, double n, int n1, vector<double> xy, vector<vector<Complex>>& funcVortex) {
	
	
	/*parallel_for_each(e,
		[=](index<2> idx) restrict(cpu, amp) {
			result[idx] = f[idx] * iexp(n * ((idx[1] < (n1 / 2)) ? atan2(-x[idx[1]], x[idx[0]]) : (atan2(-x[idx[1]], x[idx[0]]) + 2 * Pi)));
		});
	result.synchronize();*/
}

double funcGauss2D(double a, int n1, vector<double>& xy, vector<vector<double>>& func, double sigma) {
	double h = 2 * a / n1;
	double xyValue = -a;
	vector<double> func1D;

	for (int i = 0; i < n1; i++) {
		xy.push_back(xyValue);
		func1D.push_back((exp(-(xyValue * xyValue) / (2 * sigma * sigma))));
		xyValue += h;
	}

	for (int i = 0; i < n1; i++) {
		func.push_back(vector<double>());
		for (int j = 0; j < n1; j++) {
			func.at(i).push_back(func1D.at(i) * func1D.at(j));
		}
	}
	return h;
}

void collins2D(double a, double b, double A, double B, double C, double D, int n1, int n2, double n, double wavelength) {
	vector<double> xy;
	vector<vector<double>> func;

	double h;
	cout << "Входная функция:" << "\nГауссов пучок (1): ";
	string select;
	cin >> select;
	if (select == "1") {
		cout << "Введите параметр сигма:" << "\nsigma = ";
		double sigma;
		cin >> sigma;
		h = funcGauss2D(a, n1, xy, func, sigma);
	}

	vector<vector<Complex>> complex;

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

			double A, B, C, D;
			cout << "Введите коэффициенты ABCD-матрицы: " << endl;
			cin >> A >> B >> C >> D;
			if ((A * D - B * C) != 1) error("Определитель ABCD-матрицы должен быть равен 1");

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

			collins2D(a, b, A, B, C, D, n1, n2, n, wavelength);

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
