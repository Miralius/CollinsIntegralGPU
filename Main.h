#ifndef MAIN_H
#define MAIN_H

#include "field.h"
#include "IO.h"

//Функция проверки ввода и сброса потока ввода в случае неверного ввода
inline void wrongInput() {
	std::cout << "Неверный ввод! Введите ещё раз!" << std::endl;
	std::cin.clear();
	std::cin.ignore(std::cin.rdbuf()->in_avail(), '\n');
}

#endif