#include <iostream>
#include <Windows.h>

using namespace std;

inline void error(const std::string& s)
{
	throw runtime_error(s);
}

double funcGauss2D(double a, int n1, double* xy, double** func, double sigma) {
	double h = 2 * a / n1;

	return h;
}

void collins2D(double a, double b, double A, double B, double C, double D, int n1, int n2, double n, double wavelength) {
	double* xy = new double[n1];
	double** func = new double* [n1 * n1];
	for (int i = 0; i < n1; i++) {
		func[i] = new double[n1];
	}
	double h;
	cout << "������� �������:" << "\n������� ����� (1): ";
	string select;
	cin >> select;
	if (select == "1") {
		cout << "������� �������� �����:" << "\nsigma = ";
		double sigma;
		cin >> sigma;
		h = funcGauss2D(a, n1, xy, func, sigma);
	}
}

int main()
{
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

	try {
		cout << "������ ����������� ��������� ���������" << endl;
		while (1) {
			cout << "������� ������� �������������� (a � b):" << "\na = ";
			double a, b;
			cin >> a;
			if (a == 0) break;
			cout << "b = ";
			cin >> b;

			double A, B, C, D;
			cout << "������� ������������ ABCD-�������: " << endl;
			cin >> A >> B >> C >> D;
			if ((A * D - B * C) != 1) error("������������ ABCD-������� ������ ���� ����� 1");

			cout << "������� ���������� �������� �������������� (n �������� ���� � n ��������� ����):" << "\nn1 = ";
			int n1, n2;
			cin >> n1;
			cout << "n2 = ";
			cin >> n2;

			cout << "������� ����� ���������� �������� �����:" << "\nn = ";
			double n;
			cin >> n;

			cout << "������� ����� ����� �����:" << "\nwavelength = ";
			double wavelength;
			cin >> wavelength;

			collins2D(a, b, A, B, C, D, n1, n2, n, wavelength);

			cout << "���������� �������? ��� ������ ������ 0" << endl;
		}
	}
	catch (runtime_error & e) {
		cerr << "������! " << e.what() << endl;
	}
	catch (...) {
		cerr << "����������� ������!" << endl;
	}
	return 0;
}
