#ifndef IO_H
#define IO_H

#include <iostream>
#include <vector>
#include <fstream>

// �������� ����� (nameFiel � ��� �����)
template <typename T> T loadingFile(std::string nameFile) {
	std::ifstream in(nameFile, std::ios::binary | std::ios::in);
	in.unsetf(std::ios_base::skipws);
	T data{};
	if (!in.fail()) {
		if (!(in >> data))
		{
			throw std::runtime_error("���� " + nameFile + " ���� ��� �������� �������� ������!");
		}
	}
	else {
		{
			throw std::runtime_error("���� " + nameFile + " �� ������!");
		}
	}
	return data;
}

// �������� ������� �� ����� (nameFile � ��� �����)
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
			throw std::runtime_error("���� " + nameFile + " ���� ��� �������� �������� ������!");
		}
	}
	else
	{
		throw std::runtime_error("���� " + nameFile + " �� ������!");
	}
	return vectorName;
}

// ������ � ���� (data � ������ ��� ������, nameFile � ��� �����)
template <typename T> void writingFile(T& data, std::string nameFile) {
	std::ofstream output(nameFile, std::ios::binary | std::ios::trunc | std::ios::out);
	if (!output) {
		{
			throw std::runtime_error("������ � ���� " + nameFile + " ����������!");
		}
	}
	output << data;
}

#endif