#include "Main.h"

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
			while (!(cin >> limit) || !(limit > 0)) {
				wrongInput();
			}
			limits.push_back(limit);
			cout << "c = ";
			while (!(cin >> limit) || !(limit > 0)) {
				wrongInput();
			}
			limits.push_back(limit);
			cout << "d = ";
			while (!(cin >> limit) || !(limit > 0)) {
				wrongInput();
			}
			limits.push_back(limit);
			cout << "Дробное преобразование Фурье (y или n)?" << endl << "Ответ: ";
			char answer;
			while (!(cin >> answer) || !(answer == 'y' || answer == 'n')) {
				wrongInput();
			}
			vector<vector<double>> matrixABCD;
			double f, z;
			bool FrFT;
			if (answer == 'y') {
				FrFT = true;
				cout << "Введите фокусное расстояние линзы (f, мм) и расстояние до изображения (z, мм):" << endl << "f = ";
				while (!(cin >> f) || !(f > 0)) {
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
				FrFT = false;
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

			vector<vector<double>> fieldParameters;
			fieldParameters.push_back(vector<double>());

			cout << "Введите длину волны света (нм):" << endl << "wavelength = ";
			double wavelength;
			while (!(cin >> wavelength) || !(wavelength > 0)) {
				wrongInput();
			}
			fieldParameters.at(0).push_back(wavelength / 1000000);
			
			cout << "Выберите шаблон входного поля:" << endl << "Одна мода (1)" << endl << "Радиально-симметричный кластер (2): ";
			int select;
			int N;
			bool different;
			bool clusterMode;
			double parameter;
			patternField pattern;
			cin >> select;
			switch (select) {
			case 1:
				different = false;
				clusterMode = false;
				fieldParameters.at(0).push_back(1); // N = 1
				fieldParameters.at(0).push_back(0); // ro = 0
				fieldParameters.at(0).push_back(0); // angle = 0
				fieldParameters.push_back(setInputFunction(clusterMode));
				break;
			case 2:
				clusterMode = true;
				cout << "В кластере все пучки одинаковы (y или n)?" << endl << "Ответ: ";
				while (!(cin >> answer) || !(answer == 'y' || answer == 'n')) {
					wrongInput();
				}
				different = answer == 'n' ? true : false;
				cout << "Введите количество пучков:" << endl << "N = ";
				while (!(cin >> N)) {
					wrongInput();
				}
				fieldParameters.at(0).push_back(static_cast<double>(N));
				cout << "Введите радиус кластера, мм:" << endl << "ro = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.at(0).push_back(parameter);
				cout << "Последний пучок проходит через ось Ox (y или n)?" << endl << "Ответ: ";
				while (!(cin >> answer) || !(answer == 'y' || answer == 'n')) {
					wrongInput();
				}
				fieldParameters.at(0).push_back(answer == 'y' ? 1 : 0);
				if (different) {
					for (auto i = 1; i <= N; i++) {
						fieldParameters.push_back(setInputFunction(clusterMode));
					}
				}
				else {
					fieldParameters.push_back(setInputFunction(clusterMode));
				}
				break;
			default:
				error("Выбран несуществующий шаблон!");
			}
			pattern = static_cast<patternField>(select);

			if (FrFT) {
				fieldParameters.at(0).push_back(z);
			}
			else {
				cout << "Введите расстояние между входной и выходной плоскостями изображения (мм):" << endl << "z = ";
				while (!(cin >> parameter)) {
					wrongInput();
				}
				fieldParameters.at(0).push_back(parameter);
			}

			auto output = field(limits, n1, n2, crossSection::Oxy, matrixABCD, pattern);
			output.setFieldParameters(different, fieldParameters);
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

			//debug cycle
			
			/*for (auto i = 900; i <= 2000; i = (i == 1100) ? 2000 : i + 25) {
				cout << endl;
				limits.at(2) = (((i == 900) || (i == 1100)) ? 1.5 : (((i == 925) || (i == 1075)) ? 1.25 : (((i == 950) || (i == 1050)) ? 1 : (((i == 975) || (i == 1025)) ? 0.75 : 0.5))));
				limits.at(3) = (((i == 900) || (i == 1100)) ? 1.5 : (((i == 925) || (i == 1075)) ? 1.25 : (((i == 950) || (i == 1050)) ? 1 : (((i == 975) || (i == 1025)) ? 0.75 : 0.5))));
				z = static_cast<double>(i);
				fieldParameters.at(0).at(4) = z;
				matrixABCD.at(0).at(0) = cos(z / f * M_PI / 2);
				matrixABCD.at(0).at(1) = f * sin(z / f * M_PI / 2);
				matrixABCD.at(1).at(0) = -sin(z / f * M_PI / 2) / f;
				matrixABCD.at(1).at(1) = cos(z / f * M_PI / 2);
				output = field(limits, n1, n2, crossSection::Oxy, matrixABCD, pattern);
				output.setFieldParameters(different, fieldParameters);
				writingFile<BMP>(output.createBMP(absSchemeName, false), to_string(i) + ".bmp");
			}*/
			
			cout << endl << "Результаты записаны! Продолжить расчёты? Для выхода ввести 0" << endl;
		}
	}
	catch (runtime_error & e) {
		cerr << endl << "Ошибка! " << e.what() << endl;
		system("pause");
	}
	catch (...) {
		cerr << "Неизвестная ошибка!" << endl;
	}
	return 0;
}