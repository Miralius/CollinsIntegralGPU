#pragma once
#include "field.h"
class RadialCluster :
    public field
{
public:
    // ����������� ����������� �������� �������� ������ �� ���������� ���������� [begin, end] �� x � y.
    // ����� ������������� ���������� ������, ������ ��������. ������� ���������: begin � ������ ��������� (��),
    // end � ����� ��������� (��), count � ���������� ����� � ���������, modeCount � ���������� ������,
    // radius � ������ �������� (��).
    RadialCluster(double begin, double end, double radius, unsigned count, unsigned modeCount);
    // ������� �������������, �� ����� ��� ����� �������� �������.
    // ��������� � �������� ��������� ������ ������� complex<double>(double, double, unsigned).
    void setMode(std::function<std::complex<double>(double, double, unsigned)>&& mode);
    // ������� �������������, �� ����� ��� ����� �������� �������.
    // ��������� � �������� ��������� ������ ������� complex<double>(double, double, unsigned).
    // ���������� ������ ��������� ����.
    void setAndCalculateMode(std::function<std::complex<double>(double, double, unsigned)>&& mode);
private:
    const unsigned modeCount;
    const double radius;
    std::function<std::complex<double>(double, double, unsigned)> customMode; //f_n(x,y,N)
    void calculate();
};

