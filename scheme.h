#ifndef SCHEME_H
#define SCHEME_H

#include <vector>
#include <Windows.h>
#include <algorithm>
#include <iostream>
#include "IO.h"

// Класс цветовой схемы
class scheme {
private:
	std::vector<std::vector<byte>> colorScheme; // Цветовая схема
	std::string schemeName; // Имя схемы

public:
	// Конструктор по умолчанию
	scheme();
	// Конструктор цветовой схемы (schemeName — путь с именем файла до цветовой схемы)
	scheme(std::string schemeName);
	// Оператор преобразования цветовой схемы в двумерный вектор байтов
	operator std::vector<std::vector<byte>>();
};

// Класс пикселя
class pixel {
private:
	std::vector<byte> colors; // Вектор субпикселей

public:
	// Конструктор пикселя по умолчанию
	pixel();
	// Конструктор пикселя (blue — синий канал, green — зелёный, red — красный, alpha — прозрачность)
	pixel(byte blue, byte green, byte red, byte alpha);
	// Оператор преобразования пикселя в вектор субпикселей
	operator std::vector<byte>();
};

// Оператор ввода пикселей (input — поток ввода, obj — пиксель)
std::istream& operator>>(std::istream& input, pixel& obj);

#endif