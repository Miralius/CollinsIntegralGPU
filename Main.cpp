#include "Main.h"
#include <QtCore/QCoreApplication>

field getRadialSymmetricClusterMode(const field& solitone, double sigma, double radius, int number, double ksi, double eta)
{
	field superposition = solitone;
	superposition.shift(radius, 0);
	auto parameterMode = field(superposition.getY().back(), -superposition.getY().back(), static_cast<int>(superposition.getY().size()));
	parameterMode.initialTransverseVelocityAndPowerFactorExpMode(ksi, eta, sigma, radius, 0);
	superposition *= parameterMode;
	for (auto i = 1; i < number; i++) {
		field solitoneShifted = solitone;
		auto phi_n = M_PI * 2 * i / number;
		auto c_xn = radius * cos(phi_n);
		auto c_yn = radius * sin(phi_n);
		solitoneShifted.shift(c_xn, c_yn);
		parameterMode.initialTransverseVelocityAndPowerFactorExpMode(ksi, eta, sigma, c_xn, c_yn);
		solitoneShifted *= parameterMode;
		superposition += solitoneShifted;
	}
	return superposition;
}

int main(int argc, char* argv[]) {

	QCoreApplication a(argc, argv);

	return a.exec();
	/*SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		std::cout << std::endl << "Моделирование входного пучка" << std::endl;
		auto a = -10.;
		auto b = 10.;
		auto n = 500;
		auto solitone = field(a, b, n);
		auto alpha = 5.;
		auto beta = 5.;
		auto alpha0 = -5.;
		auto beta0 = -5.;
		auto sigma = 1.;
		auto ksi = 1.;
		auto eta = 1.;
		auto r = 5;
		auto number = 8;
		solitone.gaussMode(sigma, 1);
		auto aperture = field(a, b, n);
		aperture.vortexMode(0);
		solitone *= aperture;
		auto superposition = getRadialSymmetricClusterMode(solitone, sigma, r, number, ksi, eta);
		std::string absFileName = "input_abs.bmp";
		std::string argFileName = "input_phase.bmp";
		std::string absSchemeName = "fire";
		std::string argSchemeName = "grays";
		BMP test = superposition.createBMP(absSchemeName, false);
		BMP test2 = superposition.createBMP(argSchemeName, true);
		writingFile<BMP>(test, absFileName);
		writingFile<BMP>(test2, argFileName);
		auto u = 0.;
		auto wavelength = 650. / 1000000;
		auto z_begin = 50.;
		auto z_end = 1950.;
		auto z_n = 2000;
		auto f = 1000.;
		for (auto i = 750; i <= 1750; i += 250) {
			std::cout << std::endl << "Моделирование ДрПФ при z = " << i << std::endl;
			field oxy = superposition;
			oxy.ouvFractionalFourierTransform(a, b, n, wavelength, i, f);
			absFileName = "oxy_abs_" + std::to_string(i) + ".bmp";
			argFileName = "oxy_arg_" + std::to_string(i) + ".bmp";
			test = oxy.createBMP(absSchemeName, false);
			test2 = oxy.createBMP(argSchemeName, true);
			writingFile<BMP>(test, absFileName);
			writingFile<BMP>(test2, argFileName);
		}
		std::cout << std::endl << "Моделирование продольного сечения пучка" << std::endl;
		field oxz = superposition;
		oxz.ovzFractionalFourierTransform(a, b, n, z_begin, z_end, z_n, wavelength, u, f);
		absFileName = "oxz_abs.bmp";
		test = oxz.createBMP(absSchemeName, false);
		writingFile<BMP>(test, absFileName);
	}
	catch (std::runtime_error & e) {
		std::cerr << std::endl << "Ошибка! " << e.what() << std::endl;
		system("pause");
	}
	catch (...) {
		std::cerr << "Неизвестная ошибка!" << std::endl;
	}
	return 0;*/
}