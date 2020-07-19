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