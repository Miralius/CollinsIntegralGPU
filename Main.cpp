#include "Main.h"

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		std::cout << std::endl << "Моделирование входного пучка" << std::endl;
		auto a = -3.;
		auto b = 3.;
		auto n = 100;
		auto solitone = field(a, b, n);
		auto alpha = 5.;
		auto beta = 5.;
		auto alpha0 = -5.;
		auto beta0 = -5.;
		auto sigma = 1.;
		/*auto ksi = 0.;
		auto eta = 1.;
		auto r = 1.5;
		auto number = 8;*/
		solitone.airyMode(alpha, beta, alpha0, beta0);
		auto aperture = field(a, b, n);
		aperture.gaussMode(sigma, 2);
		solitone *= aperture;
		/*auto superposition = solitone;
		superposition.shift(r, 0);
		superposition.setInitialTransverseVelocityAndPowerFactor(ksi, eta, sigma, r, 0.);
		for (auto i = 1; i < number; i++) {
			auto solitoneShifted = solitone;
			auto phi_n = M_PI * 2 * i / number;
			auto c_xn = r * cos(phi_n);
			auto c_yn = r * sin(phi_n);
			solitoneShifted.shift(c_xn, c_yn);
			solitoneShifted.setInitialTransverseVelocityAndPowerFactor(ksi, eta, sigma, c_xn, c_yn);
			superposition += solitoneShifted;
		}*/
		std::string absFileName = "input_abs.bmp";
		std::string argFileName = "input_phase.bmp";
		std::string absSchemeName = "fire";
		std::string argSchemeName = "rainbow";
		BMP test = solitone.createBMP(absSchemeName, false);
		BMP test2 = solitone.createBMP(argSchemeName, true);
		writingFile<BMP>(test, absFileName);
		writingFile<BMP>(test2, argFileName);
		auto u = 0.;
		auto wavelength = 650. / 1000000;
		auto z_begin = 50.;
		auto z_end = 1950.;
		auto z_n = 2000;
		auto f = 1000.;
		for (auto i = 250; i <= 1750; i += 250) {
			std::cout << std::endl << "Моделирование ДрПФ при z = " << i << std::endl;
			field oxy = solitone;
			oxy.ouvFractionalFourierTransform(a, b, n, wavelength, i, f);
			absFileName = "oxy_abs_" + std::to_string(i) + ".bmp";
			argFileName = "oxy_arg_" + std::to_string(i) + ".bmp";
			test = oxy.createBMP(absSchemeName, false);
			test2 = oxy.createBMP(argSchemeName, true);
			writingFile<BMP>(test, absFileName);
			writingFile<BMP>(test2, argFileName);
		}
		std::cout << std::endl << "Моделирование продольного сечения пучка" << std::endl;
		field oxz = solitone;
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
	return 0;
}