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
			char answer;
			while (!(cin >> answer) || !(answer == 'y' || answer == 'n')) {
				wrongInput();
			}
			vector<vector<double>> matrixABCD;
			double f, z;
			if (answer == 'y') {
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

			vector<vector<double>> fieldParameters;
			fieldParameters.push_back(vector<double>());

			cout << "Введите длину волны света (нм):" << endl << "wavelength = ";
			double wavelength;
			while (!(cin >> wavelength)) {
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
				fieldParameters.at(0).push_back(1);
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