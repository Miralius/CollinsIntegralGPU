#ifndef MAIN_H
#define MAIN_H

#include <Windows.h>
#include "field.h"
#include "IO.h"

//Функция проверки ввода и сброса потока ввода в случае неверного ввода
inline void wrongInput() {
	std::cout << "Неверный ввод! Введите ещё раз!" << std::endl;
	std::cin.clear();
	std::cin.ignore(std::cin.rdbuf()->in_avail(), '\n');
}

field getRadialSymmetricClusterMode(const field& solitone, double sigma, double radius, int number, double ksi = 0.0, double eta = 0.0);

#endif