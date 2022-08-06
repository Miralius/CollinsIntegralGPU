#include "Main.h"

int main(int argc, char* argv[]) {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		auto a = -10.;
		auto b = 10.;
		auto r = 5.;
		auto n = 1000;
		auto N = 5;
		auto sigma = 1.;
		auto cluster = RadialCluster(a, b, r, n, N);
		
		auto doubledLaguerre = [sigma](double x, double y, unsigned N) {
			std::complex<double> sum{};
			for (unsigned p = 1; p <= N; ++p) {
				sum += 1. * std::pow(std::complex<double>(x, y), p);
			}
			return std::complex<double>(exp(-(x * x + y * y) / (sigma * sigma))) * sum;
		};

		std::cout << "Моделирование входного пучка " << std::endl;
		cluster.setAndCalculateMode(std::move(doubledLaguerre));
		std::string absSchemeName = "fire";
		std::string argSchemeName = "grays";
		auto u = 0.;
		auto wavelength = 650. / 1000000;
		auto z_begin = 0.;
		auto z_end = 2000.;
		auto z_n = 2000;
		auto f = 1000.;
		for (auto i = 0; i <= 2000; i += 50) {
			std::cout << "Моделирование ДрПФ пучка при z = " << i << std::endl;
			std::string absFileName = "Modes\\Cluster\\oxy_abs_" + std::to_string(i) + ".bmp";
			std::string argFileName = "Modes\\Cluster\\oxy_arg_" + std::to_string(i) + ".bmp";
			if (std::filesystem::exists(absFileName)) {
				std::cout << "Данный пучок уже смоделирован!" << std::endl;
				continue;
			}
			RadialCluster oxy = cluster;
			oxy.ouvFractionalFourierTransform(a, b, n, wavelength, i, f);
			BMP test = oxy.createBMP(absSchemeName, false);
			BMP test2 = oxy.createBMP(argSchemeName, true);
			writingFile<BMP>(test, absFileName);
			writingFile<BMP>(test2, argFileName);
		}
		std::cout << "Моделирование продольного сечения пучка" << std::endl;
		field oxz = cluster;
		std::string absFileName = "Modes\\Cluster\\oxz_abs.bmp";
		oxz.ovzFractionalFourierTransform(a, b, n, z_begin, z_end, z_n, wavelength, u, f);
		BMP test = oxz.createBMP(absSchemeName, false);
		writingFile<BMP>(test, absFileName);
	}
	catch (std::runtime_error & e) {
		std::cerr << std::endl << "Ошибка! " << e.what() << std::endl;
		system("pause");
	}
	catch (...) {
		std::cerr << "Неизвестная ошибка!" << std::endl;
	}
	return 0;
}