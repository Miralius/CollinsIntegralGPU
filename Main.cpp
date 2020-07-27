#include "CollinsIntegralGPU.h"

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		cout << "Расчёт двумерного интеграла Коллинза…" << endl;

		while (1) {
			vector<double> limits;
			double limit;
			cout << "Введите пределы интегрирования (a, b и c, d, мм):" << "\na = ";
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
				cout << "Введите фокусное расстояние линзы (f, мм) и расстояние до изображения (z, мм):" << "\nf = ";
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
			
			cout << "Входная функция:" << "\nМода Гаусса (1)\nМода Гаусса-Эрмита (2)\nМода Гаусса-Лагерра (3)\nМода Гаусса для кластерного режима (4): ";
			int select;
			double parameter;
			inputField selectedInputField;
			cin >> select;
			selectedInputField = static_cast<inputField>(select);
			switch (select) {
			case 1:
				cout << "Введите ширину пучка (мм):" << "\nsigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите топологический заряд:" << "\nm = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter);
				break;
			case 2:
				cout << "Введите ширину пучка (мм):" << "\nsigma = ";
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
				cout << "Введите ширину пучка (мм):" << "\nsigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите порядок m:" << "\nm = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите порядок n:" << "\nn = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter);
				break;
			case 4:
				cout << "Введите ширину пучка (мм):" << "\nsigma = ";
				while (!(cin >> parameter) || !(parameter > 0)) wrongInput();
				fieldParameters.push_back(parameter);
				cout << "Введите топологический заряд:" << "\nm = ";
				while (!(cin >> parameter)) wrongInput();
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
				cout << "«Кластерный режим» (y)?" << "\nОтвет: ";
				while (!(cin >> mode) || !(mode == 'y')) {
					wrongInput();
				}
			}
				
			if (mode == 'y') {
				int N;
				cout << "Введите количество пучков:" << "\nN = ";
				while (!(cin >> N)) wrongInput();
				fieldParameters.push_back(static_cast<double>(N));
				cout << "Введите радиус k*sigma:" << "\nk = ";
				while (!(cin >> parameter)) wrongInput();
				fieldParameters.push_back(parameter * fieldParameters.at(1));
			}

			auto output = field(limits, n1, n2, crossSection::Oxy, matrixABCD, selectedInputField, fieldParameters);
			if (mode == 'y') {
				output.setClusterMode(true);
			}
			string outputAbs;
			string outputArg;
			cout << "Введите название файла для амплитуды: ";
			cin >> outputAbs;
			cout << "Введите название файла для фазы: ";
			cin >> outputArg;
			writingFile<BMP>(output.createBMP("fire", false), outputAbs);
			writingFile<BMP>(output.createBMP("black_white", true), outputArg);
			
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