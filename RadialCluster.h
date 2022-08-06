#pragma once
#include "field.h"
class RadialCluster :
    public field
{
public:
    // Конструктор радиального кластера световых пучков по параметрам интервалов [begin, end] по x и y.
    // Также устанавливает количество пучков, радиус кластера. Входные параметры: begin — начало интервала (мм),
    // end — конец интервала (мм), count — количество точек в интервале, modeCount — количество пучков,
    // radius — радиус кластера (мм).
    RadialCluster(double begin, double end, double radius, unsigned count, unsigned modeCount);
    // Функция устанавливает, из каких мод будет состоять кластер.
    // Принимает в качестве аргумента лямбда функцию complex<double>(double, double, unsigned).
    void setMode(std::function<std::complex<double>(double, double, unsigned)>&& mode);
    // Функция устанавливает, из каких мод будет состоять кластер.
    // Принимает в качестве аргумента лямбда функцию complex<double>(double, double, unsigned).
    // Производит расчёт светового поля.
    void setAndCalculateMode(std::function<std::complex<double>(double, double, unsigned)>&& mode);
private:
    const unsigned modeCount;
    const double radius;
    std::function<std::complex<double>(double, double, unsigned)> customMode; //f_n(x,y,N)
    void calculate();
};

