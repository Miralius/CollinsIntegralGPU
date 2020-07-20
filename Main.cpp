#include "CollinsIntegralGPU.h"

extern int test();

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	test();

	try {
		cout << "Расчёт двумерного интеграла Коллинза…" << endl;

		while (1) {
			vector<double> limits;
			double limit;
			cout << "Введите пределы интегрирования (a, b и c, d):" << "\na = ";
			while (!(cin >> limit)) {
				wrongInput();
			}
			limits.push_back(limit);
			if (limits.at(0) == 0) {
				break;
			}
			cout << "b = ";
			while (!(cin >> limit)) {
				wrongInput();
			}
			limits.push_back(limit);
			cout << "c = ";
			while (!(cin >> limit)) {
				wrongInput();
			}
			limits.push_back(limit);
			cout << "d = ";
			while (!(cin >> limit)) {
				wrongInput();
			}
			limits.push_back(limit);
			cout << "Дробное преобразование Фурье (y или n)?" << "\nОтвет: ";
			char fourier;
			while (!(cin >> fourier) || !(fourier == 'y' || fourier == 'n')) {
				wrongInput();
			}
			vector<vector<double>> matrixABCD;
			if (fourier == 'y') {
				double f, z;
				cout << "Введите фокусное расстояние линзы (f, мм) и расстояние до изображения:" << "\nf = ";
				while (!(cin >> f)) {
					wrongInput();
				}
				cout << "z = ";
				while (!(cin >> z)) {
					wrongInput();
				}
				matrixABCD.push_back(vector<double>());
				matrixABCD.at(0).push_back(cos(z / f * M_PI / 2));
				matrixABCD.at(0).push_back(f * sin(z / f * M_PI / 2));
				matrixABCD.push_back(vector<double>());
				matrixABCD.at(1).push_back(-sin(z / f * M_PI / 2) / f);
				matrixABCD.at(1).push_back(cos(z / f * M_PI / 2));
			}
			else {
				cout << "Введите коэффициенты ABCD-матрицы: " << endl;
				double coefficient;
				for (auto i = 0; i < 2; i++) {
					matrixABCD.push_back(vector<double>());
					for (auto j = 0; j < 2; j++) {
						while (!(cin >> coefficient)) {
							wrongInput();
						}
						matrixABCD.back().push_back(coefficient);
					}
				}
				if (std::abs((matrixABCD.at(0).at(0) * matrixABCD.at(1).at(1) - matrixABCD.at(0).at(1) * matrixABCD.at(1).at(0)) - 1) > FLT_EPSILON) {
					error("Определитель ABCD-матрицы должен быть равен 1");
				}
			}
			
			cout << "Введите количество отсчётов интегрирования (n входного поля и n выходного поля):" << "\nn1 = ";
			int n1, n2;
			while (!(cin >> n1) || !(n1 > 0)) {
				wrongInput();
			}
			cout << "n2 = ";
			while (!(cin >> n2) || !(n2 > 0)) {
				wrongInput();
			}

			vector<double> fieldParameters;

			cout << "Введите длину волны света (нм):" << "\nwavelength = ";
			double wavelength;
			while (!(cin >> wavelength)) {
				wrongInput();
			}
			fieldParameters.push_back(wavelength / 1000000);
			
			cout << "Входная функция:" << "\nМода Гаусса (1)\nМода Гаусса-Эрмита (2)\nМода Гаусса-Лагерра (n = 0, m) (3): ";
			int select;
			double parameter;
			inputField selectedInputField;
			cin >> select;
			selectedInputField = static_cast<inputField>(select);
			switch (select) {
			case 1:
				cout << "Введите параметр сигма:" << "\nsigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите топологический заряд:" << "\nm = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter);
				break;
			case 2:
				cout << "Введите параметр сигма:" << "\nsigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите порядок m:" << "\nm = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите порядок n:" << "\nn = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter);
				break;
			case 3:
				cout << "Введите параметр сигма:" << "\nsigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите порядок m:" << "\nm = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите порядок n:" << "\nn = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter);
				break;
			default:
				error("Не выбрана исходная функция!");
			}

			auto output = field(limits, n1, n2, crossSection::Oxy, matrixABCD, selectedInputField, fieldParameters);
			writingFile<BMP>(output.createBMP("fire", false), "absOutput.bmp");
			writingFile<BMP>(output.createBMP("black_white", true), "argOutput.bmp");

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