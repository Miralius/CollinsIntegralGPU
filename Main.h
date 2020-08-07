#ifndef MAIN_H
#define MAIN_H

#include "field.h"

using namespace std;

inline void error(const string& s) {
	throw runtime_error(s);
}

inline void wrongInput() {
	cout << "Неверный ввод! Введите ещё раз!" << endl;
	cin.clear();
	cin.ignore(cin.rdbuf()->in_avail(), '\n');
}

inline vector<double> setInputFunction(bool clusterMode) {
	cout << "Входная функция:" << endl << "Мода Гаусса (1)" << endl << "Мода Гаусса-Эрмита (2)" << endl << "Мода Гаусса-Лагерра (3): ";
	int select;
	double parameter;
	vector<double> modeParameters;
	cin >> select;
	modeParameters.push_back(static_cast<double>(select));
	switch (select) {
	case 1:
		cout << "Введите ширину пучка (мм):" << endl << "beam waist = ";
		while (!(cin >> parameter) || !(parameter > 0)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		cout << "Введите значение коэффициента перед шириной пучка (обычно 2):" << endl << "beam waist coefficient = ";
		while (!(cin >> parameter) || !(parameter > 0)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		if (clusterMode) {
			cout << "Начальная фаза равна:" << endl << "exp(imф) (0)" << endl << "exp(iфj) (1): ";
			while (!(cin >> parameter)) {
				wrongInput();
			}
		}
		else {
			parameter = 0;
		}
		modeParameters.push_back(parameter);
		cout << "Введите топологический заряд:" << endl << "m = ";
		while (!(cin >> parameter)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		break;
	case 2:
		cout << "Введите ширину пучка (мм):" << endl << "beam waist = ";
		while (!(cin >> parameter) || !(parameter > 0)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		cout << "Введите значение коэффициента перед шириной пучка (обычно 2):" << endl << "beam waist coefficient = ";
		while (!(cin >> parameter) || !(parameter > 0)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		modeParameters.push_back(0);
		cout << "Введите порядок m:" << endl << "m = ";
		while (!(cin >> parameter)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		cout << "Введите порядок n:" << endl << "n = ";
		while (!(cin >> parameter)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		break;
	case 3:
		cout << "Введите ширину пучка (мм):" << endl << "beam waist = ";
		while (!(cin >> parameter) || !(parameter > 0)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		cout << "Введите значение коэффициента перед шириной пучка (обычно 2):" << endl << "beam waist coefficient = ";
		while (!(cin >> parameter) || !(parameter > 0)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		modeParameters.push_back(0);
		cout << "Введите порядок m:" << endl << "m = ";
		while (!(cin >> parameter)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		cout << "Введите порядок n:" << endl << "n = ";
		while (!(cin >> parameter)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		break;
	default:
		error("Не выбрана исходная функция!");
	}
	if (clusterMode) {
		while (modeParameters.size() < 6) {
			modeParameters.push_back(0);
		}
		cout << "Введите несущую (начальную поперечную) скорость пучка (мм/c):" << endl << "ksi = ";
		while (!(cin >> parameter)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
		cout << "Введите коэффициент мощности пучка:" << endl << "eta = ";
		while (!(cin >> parameter)) {
			wrongInput();
		}
		modeParameters.push_back(parameter);
	}
	return modeParameters;
}

template <typename T> T loadingFile(string nameFile) {
	ifstream in(nameFile, ios::binary | ios::in);
	in.unsetf(ios_base::skipws);
	T data;
	if (!in.fail()) {
		if (!(in >> data)) error("Файл " + nameFile + " пуст или содержит неверные данные!");
	}
	else {
		error("Файл " + nameFile + " не найден!");
	}
	return data;
}

template <typename T> vector<T> loadingData(string nameFile) {
	ifstream in(nameFile);
	vector<T> vectorName;
	if (!in.fail()) {
		T buffer;
		while (in >> buffer) {
			if (in.eof()) break;
			vectorName.push_back(buffer);
		}
		in.close();
		if (vectorName.size() == 0) error("Файл " + nameFile + " пуст или содержит неверные данные!");
	}
	else error("Файл " + nameFile + " не найден!");
	return vectorName;
}

template <typename T> void writingFile(T& data, string nameFile) {
	ofstream output(nameFile, ios::binary | ios::trunc | ios::out);
	if (!output) {
		error("Запись в файл " + nameFile + " невозможна!");
	}
	output << data;
}

#endif