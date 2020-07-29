#include "CollinsIntegralGPU.h"

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		cout << "Расчёт двумерного интеграла Коллинза…" << endl;

		while (1) {
			vector<double> limits;
			double limit;
			cout << "Введите пределы интегрирования (a, b и c, d, мм):" << endl << "a = ";
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
			cout << "Дробное преобразование Фурье (y или n)?" << endl << "Ответ: ";
			char fourier;
			while (!(cin >> fourier) || !(fourier == 'y' || fourier == 'n')) {
				wrongInput();
			}
			vector<vector<double>> matrixABCD;
			double f, z;
			if (fourier == 'y') {
				cout << "Введите фокусное расстояние линзы (f, мм) и расстояние до изображения (z, мм):" << endl << "f = ";
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
			
			cout << "Введите количество отсчётов интегрирования (n входного поля и n выходного поля):" << endl << "n1 = ";
			int n1, n2;
			while (!(cin >> n1) || !(n1 > 0)) {
				wrongInput();
			}
			cout << "n2 = ";
			while (!(cin >> n2) || !(n2 > 0)) {
				wrongInput();
			}

			vector<double> fieldParameters;

			cout << "Введите длину волны света (нм):" << endl << "wavelength = ";
			double wavelength;
			while (!(cin >> wavelength)) {
				wrongInput();
			}
			fieldParameters.push_back(wavelength / 1000000);
			
			cout << "Входная функция:" << endl << "Мода Гаусса (1)" << endl << "Мода Гаусса-Эрмита (2)" << endl << "Мода Гаусса-Лагерра (3)" << endl << "Мода Гаусса для кластерного режима (4): ";
			int select;
			double parameter;
			inputField selectedInputField;
			cin >> select;
			selectedInputField = static_cast<inputField>(select);
			switch (select) {
			case 1:
				cout << "Введите ширину пучка (мм):" << endl << "sigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				cout << "Введите топологический заряд:" << endl << "m = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				break;
			case 2:
				cout << "Введите ширину пучка (мм):" << endl << "sigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				cout << "Введите порядок m:" << endl << "m = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				cout << "Введите порядок n:" << endl << "n = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				break;
			case 3:
				cout << "Введите ширину пучка (мм):" << endl << "sigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				cout << "Введите порядок m:" << endl << "m = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				cout << "Введите порядок n:" << endl << "n = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				break;
			case 4:
				cout << "Введите ширину пучка (мм):" << endl << "sigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				cout << "Введите топологический заряд:" << endl << "m = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter);
				break;
			default:
				error("Не выбрана исходная функция!");
			}

			char mode;
			if (select == 4) {
				mode = 'y';
			}
			else {
				cout << "«Кластерный режим» (y)?" << endl << "Ответ: ";
				while (!(cin >> mode)) {
					wrongInput();
				}
			}
				
			if (mode == 'y') {
				if (fourier == 'y') {
					fieldParameters.push_back(z);
				}
				else {
					cout << "Введите расстояние от начальной до выходной плоскости изображения (z, мм):" << endl << "z = ";
					while (!(cin >> parameter)) {
						wrongInput();
					}
					fieldParameters.push_back(parameter);
				}
				int N;
				cout << "Введите количество пучков:" << endl << "N = ";
				while (!(cin >> N)) {
					wrongInput();
				}
				fieldParameters.push_back(static_cast<double>(N));
				cout << "Введите радиус k*sigma:" << endl << "k = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.push_back(parameter * fieldParameters.at(1));
			}

			auto output = field(limits, n1, n2, crossSection::Oxy, matrixABCD, selectedInputField, fieldParameters);
			if (mode == 'y') {
				output.setClusterMode(true);
			}
			string absFileName;
			string argFileName;
			string absSchemeName;
			string argSchemeName;
			cout << "Введите название файла для амплитуды: ";
			cin >> absFileName;
			cout << "Введите название файла для фазы: ";
			cin >> argFileName;
			cout << "Введите название цветовой схемы для амплитуды: ";
			cin >> absSchemeName;
			cout << "Введите название цветовой схемы для фазы: ";
			cin >> argSchemeName;
			writingFile<BMP>(output.createBMP(absSchemeName, false), absFileName);
			writingFile<BMP>(output.createBMP(argSchemeName, true), argFileName);
			
			cout << endl << "Результаты записаны! Продолжить расчёты? Для выхода ввести 0" << endl;
		}
	}
	catch (runtime_error & e) {
		cerr << endl << "Ошибка! " << e.what() << endl;
		main();
	}
	catch (...) {
		cerr << "Неизвестная ошибка!" << endl;
	}
	return 0;
}