#ifndef FIELD_H
#define FIELD_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iomanip>
#include <functional>
#include "scheme.h"
#include "BMP.h"

//Функция вычисления интеграла Коллинза с помощью технологии CUDA. Входные параметры: input — входное поле, x1 — вектор x, x2 — вектор y, x3 — вектор u или z, x4 — вектор v, parameters — массив параметров преобразования, dimension — массив размерностей, transformType — вид преобразования: 0 — дробное преобразование Фурье, 1 — преобразование Френеля, 2 — другое (с заданной определенной ABCD матрицей). Выходной параметр: output — выходное поле (результат вычисления интеграла Коллинза).
extern std::vector<std::vector<std::complex<double>>> calculateCollinsCUDA(const std::vector<std::vector<std::complex<double>>>& input, const std::vector<double>& x1, const std::vector<double>& x2, const std::vector<double>& x3, const std::vector<double>& x4, const std::vector<double>& parameters, const std::vector<int>& dimension, int transformType);

//Класс светового поля
class field {
private:
	std::vector<double> x; // Вектор x (до преобразования) или u (после преобразования), или z (в случае продольного сечения преобразования)
	std::vector<double> y; // Вектор y (до преобразования) или v (после преобразования)
	std::vector<std::vector<std::complex<double>>> calculatedField; // Матрица светового поля
	// Функция построения интервала (аналог linspace в Numpy). Входные параметры: begin — начало интервала (мм), end — конец интервала (мм), count — количество точек в интервале. Возвращается вектор интервала с count точками.
	static std::vector<double> calcPoints(double begin, double end, int count);
	// createField создает световое поле по заданному правилу mode (по сути в аргументе передается лямбда-функция f(x, y)). Входной параметр — лямбда-функция входного поля f(x,y). Результат записывается в поле calculatedField.
	void createField(const std::function<std::complex<double>(double, double)>& mode);
	// Статическая функция нахождения максимального по модулю значения в двумерной матрице вещественных чисел
	static double maximum(const std::vector<std::vector<double>>& field);
	// Статическая функция нахождения максимального по модулю значения в векторе комплексных чисел
	static double maximum(std::vector<std::complex<double>>& column);
	// Транспонирование матрицы. Возвращается транспонированная двумерная комплексная матрица. Входной параметр matrix — двумерная комплексная матрица, которую необходимо транспонировать.
	std::vector<std::vector<std::complex<double>>> transpose(std::vector<std::vector<std::complex<double>>>& matrix);
public:
	// Конструктор по умолчанию field
	field();
	// Конструктор, создающий экземпляр field по параметрам интервалов [a, b] по x и y. Входные параметры: a — начало интервала (мм), b — конец интервала (мм), n — количество точек в интервале.
	field(double a, double b, int n);
	// Конструктор, создающий экземпляр field, записывая в поля готовые матрицы светового поля и векторов x, y. Входные параметры: field — матрица светового поля, x — вектор x, y — вектор y.
	field(const std::vector<std::vector<std::complex<double>>>& field, const std::vector<double>& x, const std::vector<double>& y);
	// Геттер вектора x
	std::vector<double> getX();
	// Геттер вектора y
	std::vector<double> getY();
	// Геттер матрицы светового поля
	std::vector<std::vector<std::complex<double>>> getCalculatedField();
	// Сеттер матрицы светового поля
	void setCalculatedField(std::vector<std::vector<std::complex<double>>>&& calculatedField);
	// Транспонирование поля.
	void transpose();
	// Нормализация светового поля (максимальная интенсивность в каждом столбце равна 1)
	void normalize();
	// Сдвиг светового поля по координатам x и y. Входные параметры: x — сдвиг по оси x в миллиметрах (положительное значение — сдвиг вправо, отрицательное — влево), y — сдвиг по оси y в миллиметрах (положительное значение — сдвиг вверх, отрицательное — вниз).
	void shift(double x, double y);
	// Поворот светового поля. Входной параметр angle (в радианах) — поворот по часовой стрелке, если angle > 0, против часовой стрелке если angle < 0. 
	void rotate(double angle);
	// Расчёт интеграла Коллинза в случае его вырождения. Входные параметры: u — вектор u, v — вектор v, wavelength — длина волны в миллиметрах, matrixCD — матрица, состоящая только из элементов C и D.
	void collinsSingular(const std::vector<double>& u, const std::vector<double>& v, double wavelength, const std::vector<double>& matrixCD);
	// Поиск сингулярных значений в заданном промежутке z, возвращает true в случае наличия таковых, false в противном случае. Входные параметры: z — вектор z, transformType — тип преобразования (0 — дробное преобразование Фурье, 1 — преобразование Френеля, 2 — произвольное преобразование), parameter — необязательный параметр преобразования (фокусное расстояние f в мм для дробного преобразования Фурье или элемент B ABCD-матрицы в случае произвольного преобразования).
	static bool hasSingularities(const std::vector<double>& z, int transformType, double parameter = 0);
	// Расчёт интеграла Коллинза в случае поперечного сечения дробного преобразования Фурье. Входные параметры: a — начало интервала по u, v в мм, b — конец интервала по u, v в мм, n — количество точек в этих интервалах, wavelength — длина волны в мм, z — расстояние до изображения в мм, f — фокусное расстояние в мм.
	void ouvFractionalFourierTransform(double a, double b, int n, double wavelength, double z, double f);
	// Расчёт интеграла Коллинза в случае поперечного сечения преобразования Френеля. Входные параметры: a — начало интервала по u, v в мм, b — конец интервала по u, v в мм, n — количество точек в этих интервалах, wavelength — длина волны в мм, z — расстояние до изображения в мм.
	void ouvFresnelTransform(double a, double b, int n, double wavelength, double z);
	// Расчёт интеграла Коллинза в случае поперечного сечения произвольного преобразования с заданной ABCD-матрицей. Входные параметры: a — начало интервала по u, v в мм, b — конец интервала по u, v в мм, n — количество точек в этих интервалах, wavelength — длина волны в мм, z — расстояние до изображения в мм, matrixABCD — ABCD-матрица.
	void ouvTransform(double a, double b, int n, double wavelength, double z, const std::vector<double>& matrixABCD);
	// Расчёт интеграла Коллинза в случае продольного сечения дробного преобразования Фурье. Входные параметры: a — начало интервала по v в мм, b — конец интервала по v в мм, n — количество точек в интервале по v, a_z — начало интервала по z в мм, b_z — конец интервала по z в мм, n_z — количество точек в интервале по z, wavelength — длина волны в мм, u — координата сечения по u в мм, f — фокусное расстояние в мм.
	void ovzFractionalFourierTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u, double f);
	// Расчёт интеграла Коллинза в случае продольного сечения преобразования Френеля. Входные параметры: a — начало интервала по v в мм, b — конец интервала по v в мм, n — количество точек в интервале по v, a_z — начало интервала по z в мм, b_z — конец интервала по z в мм, n_z — количество точек в интервале по z, wavelength — длина волны в мм, u — координата сечения по u в мм.
	void ovzFresnelTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u);
	// Расчёт интеграла Коллинза в случае продольного сечения произвольного преобразования с заданной ABCD-матрицей. Входные параметры: a — начало интервала по v в мм, b — конец интервала по v в мм, n — количество точек в интервале по v, a_z — начало интервала по z в мм, b_z — конец интервала по z в мм, n_z — количество точек в интервале по z, wavelength — длина волны в мм, u — координата сечения по u в мм, matrixABCD — ABCD-матрица.
	void ovzTransform(double a, double b, int n, double a_z, double b_z, int n_z, double wavelength, double u, const std::vector<double>& matrixABCD);
	// Статическая функция расчёта амплитуды светового поля. Входный параметр: field — световое поле. Возвращается двумерный вектор амплитуд в каждой точке светового поля.
	static std::vector<std::vector<double>> abs(const field& field);
	// Статическая функция расчёта фазы светового поля. Входный параметр: field — световое поле. Возвращается двумерный вектор фаз в каждой точке светового поля.
	static std::vector<std::vector<double>> arg(const field& field);
	// Создает экземпляр изображения в BMP формате. Входные параметры: schemeName — название световой схемы (помещаются в папку schemes в корне программы), phase — false — строит изображение амплитуды, true — фазы.
	BMP createBMP(std::string schemeName, bool phase);
	// Задает световое поле с функций Гаусса. Входные параметры: sigma — радиус пучка, положительное значение в миллиметрах, factorSigma — коэффициент перед радиусом пучка.
	void gaussMode(double sigma, double factorSigma);
	// Задает световое поле с функцией Гаусса-Эрмита. Входные параметры: m — порядок m, n — порядок n, sigma — радиус пучка, положительное значение в миллиметрах.
	void gaussHermiteMode(int m, int n, double sigma);
	// Задает световое поле с функцией Гаусса-Лагерра. Входные параметры: m — порядок m, n — порядок n, sigma — радиус пучка, положительное значение в миллиметрах.
	void gaussLaguerreMode(int m, int n, double sigma);
	// Задает световое поле с суммой двух мод Лагерра-Гаусса. Входные параметры: sigma — радиус пучка, положительное значение в миллиметрах.
	void doubledLaguerreMode(double sigma);
	// Задает световое поле с функцией Эйри-Гаусса. Входные параметры: alpha, beta — константы, отвечающие за интенсивность пучка. alpha0 — смещение пучка по оси x (положительное значение — вправо, отрицательное — влево), beta0 — смещение пучка по оси y (положительное значение — вверх, отрицательное — вниз).
	void airyMode(double alpha, double beta, double alpha0, double beta0);
	// Задает вихрь. Входной параметр m — топологический заряд вихря.
	void vortexMode(double m);
	// Задает начальные поперечные скорости и мощности пучкам в радиально-симметричном кластере. Входные параметры: ksi — начальная поперечная скорость пучка в мм/с, eta — коэффициент мощности пучка в мВт, sigma — радиус пучка в мм, x_shift — смещение пучка относительно центра по оси x в мм (положительное значение — вправо, отрицательное — влево), y_shift — смещение пучка относительно центра по оси y в мм (положительное значение — вверх, отрицательное — вниз).
	void initialTransverseVelocityAndPowerFactorExpMode(double ksi, double eta, double sigma, double x_shift, double y_shift);
};

// Оператор суперпозиции двух световых полей. Первое поле — left, второе поле — right.
field operator+(const field& left, const field& right);
// Суперпозиция двух световых полей с присвоением. Световое поле до знака += — left, после — right.
field operator+=(field& left, const field& right);
// Поэлементное перемножение двух световых полей. Первое поле — left, второе поле — right.
field operator*(const field& left, const field& right);
// Поэлементное умножение двух световых полей с присвоением. Световое поле до знака *= — left, после — right.
field operator*=(field& left, const field& right);
#endif