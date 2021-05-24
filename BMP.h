#ifndef BMP_H
#define BMP_H

#include <vector>
#include <Windows.h>
#include <algorithm>
#include <iostream>

// Класс изображения BMP
class BMP {
private:
	static auto const countRGBChannel = 4; // Количество цветовых каналов
	static auto const BMPFILEHEADERsize = 14; // Размер заголовка BMP файла
	static auto const BMPINFOHEADERsize = 124; // Размер заголовка информации об изображении
	static auto const COLORPROFILEsize = 12; // Размер цветового профиля
	std::vector<std::vector<std::vector<byte>>> pixels; // Трёхмерная матрица изображения (первая размерность — строка, 2-ая — столбец, 3-я — вектор пикселя, значения каждого — интенсивность цвета)
	std::vector<std::vector<int>> bmpFileHeader; // Этот вектор заменяет структуру BITMAPFILEHEADER. Первый элемент каждого вложенного вектора — значение, а второй — размер типа значения.
	std::vector<std::vector<int>> bmpInfoHeader; // Этот вектор заменяет структуру BITMAPINFOHEADER. Первый элемент каждого вложенного вектора — значение, а второй — размер типа значения.
	std::vector<std::vector<int>> colorProfile; // Вектор цветового профиля
	// Функция инициализации заголовков (int width — ширина изображения в пикселях, int height — высота)
	void initHeaders(int width, int height);

public:
	// Конструктор BMP по умолчанию
	BMP();
	// Конструктор BMP (const vector<vector<vector<byte>>>& picture — трехмёрная матрица изображения (первая размерность — строка, 2-ая — столбец, 3-я — вектор пикселя, значения каждого — интенсивность цвета))
	BMP(const std::vector<std::vector<std::vector<byte>>>& picture);
	// Конструктор копирования BMP (const BMP& obj — изображение, которое необходимо скопировать)
	BMP(const BMP& obj);
	// Статический метод преобразования значения в бинарный вид (const vector<int>& number — значение (пара значение-размер значения))
	static std::vector<byte> toBinary(const std::vector<int>& number);
	// Статический метод преобразования байтов в значение (const vector<byte>& binary — набор байтов)
	static int toNumber(const std::vector<byte>& binary);
	// Оператор присваивания (const BMP& obj — объект, который необходимо присвоить)
	BMP& operator=(const BMP& obj);
	// Оператор преобразования изображения в вектор байтов
	operator std::vector<byte>();
};

// Оператор ввода изображения (input — поток ввода, bmp — ссылка на объект, в который изображение должно быть записано)
std::istream& operator>>(std::istream& input, BMP& bmp);
// Оператор вывода (output — поток вывода, bmp — изображение)
std::ostream& operator<<(std::ostream& output, const BMP& bmp);

#endif