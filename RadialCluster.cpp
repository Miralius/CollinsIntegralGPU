#include "RadialCluster.h"

RadialCluster::RadialCluster(double begin, double end, double radius, unsigned count, unsigned modeCount)
	: field(begin, end, count)
	, radius(radius)
	, modeCount(modeCount)
{
}

void RadialCluster::setMode(std::function<std::complex<double>(double, double, unsigned)>&& mode)
{
	customMode = std::move(mode);
}

void RadialCluster::setAndCalculateMode(std::function<std::complex<double>(double, double, unsigned)>&& mode)
{
	setMode(std::move(mode));
	calculate();
}

void RadialCluster::calculate()
{
	std::vector<std::vector<std::complex<double>>> clusterField;
	auto x = getX();
	auto y = getY();
	clusterField.reserve(y.size());
	auto vectorRow = std::vector<std::complex<double>>();
	vectorRow.resize(x.size());
	for (auto row : y) {
		clusterField.emplace_back(vectorRow);
	}
	for (unsigned N = 1; N <= modeCount; ++N) {
		auto phi_n = 2 * M_PI * N / modeCount;
		for (size_t row = 0; row < clusterField.size(); ++row) {
			for (size_t column = 0; column < clusterField.at(0).size(); ++column) {
				auto x_arg = x.at(column) - radius * cos(phi_n);
				auto y_arg = y.at(row) - radius * sin(phi_n);
				clusterField.at(row).at(column) += customMode(x_arg, y_arg, N);
			}
		}
	}
	setCalculatedField(std::move(clusterField));
}
