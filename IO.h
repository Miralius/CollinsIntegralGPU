#ifndef IO_H
#define IO_H

#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>

// Загрузка файла (nameFiel — имя файла)
template <typename T> T loadingFile(std::string nameFile) {
	std::ifstream in(nameFile, std::ios::binary | std::ios::in);
	in.unsetf(std::ios_base::skipws);
	T data{};
	if (!in.fail()) {
		if (!(in >> data))
		{
			throw std::runtime_error("Файл " + nameFile + " пуст или содержит неверные данные!");
		}
	}
	else {
		{
			throw std::runtime_error("Файл " + nameFile + " не найден!");
		}
	}
	return data;
}

// Загрузка вектора из файла (nameFile — имя файла)
template <typename T> std::vector<T> loadingData(std::string nameFile) {
	std::ifstream in(nameFile);
	std:: vector<T> vectorName;
	if (!in.fail()) {
		T buffer{};
		while (in >> buffer) {
			if (in.eof()) break;
			vectorName.emplace_back(buffer);
		}
		in.close();
		if (vectorName.size() == 0)
		{
			throw std::runtime_error("Файл " + nameFile + " пуст или содержит неверные данные!");
		}
	}
	else
	{
		throw std::runtime_error("Файл " + nameFile + " не найден!");
	}
	return vectorName;
}

// Запись в файл (data — данные для записи, nameFile — имя файла)
template <typename T> void writingFile(T& data, std::string nameFile) {
	if (nameFile.find('\\') != -1) {
		std::filesystem::path pathFile = nameFile;
		std::filesystem::create_directories(pathFile.remove_filename());
	}
	std::ofstream output(nameFile, std::ios::binary | std::ios::trunc | std::ios::out);
	if (!output) {
		{
			throw std::runtime_error("Запись в файл " + nameFile + " невозможна!");
		}
	}
	output << data;
}

#endif