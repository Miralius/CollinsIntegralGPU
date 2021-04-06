#include "Main.h"

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		auto a = 3.;
		auto b = 3.;
		auto n = 500;
		auto solitone = field(a, b, n);
		auto alpha = 1.;
		auto beta = 1.;
		auto alpha0 = 0.;
		auto beta0 = 0.;
		auto sigma = 0.5;
		auto ksi = 1.;
		auto eta = 1.;
		auto r = 1.5;
		auto number = 8;
		solitone.airyMode(alpha, beta, alpha0, beta0, sigma);
		auto aperture = field(a, b, n);
		aperture.gaussMode(sigma, 0.);
		solitone *= aperture;
		auto superposition = solitone;
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
		}
		string absFileName = "input_abs.bmp";
		string argFileName = "input_phase.bmp";
		string absSchemeName = "fire";
		string argSchemeName = "rainbow";
		writingFile<BMP>(superposition.createBMP(absSchemeName, false), absFileName);
		writingFile<BMP>(superposition.createBMP(argSchemeName, true), argFileName);
		auto u = 0.;
		auto wavelength = 650. / 1000000;
		auto z_begin = 100.;
		auto z_end = 1900.;
		auto z_n = 2000.;
		auto f = 1000.;
		for (auto i = 250; i <= 1750; i += 250) {
			auto oxy = superposition;
			oxy.transform(a, b, n, wavelength, transformType::fractionalFourier, i, f);
			absFileName = "oxy_abs_" + to_string(i) + ".bmp";
			argFileName = "oxy_arg_" + to_string(i) + ".bmp";
			writingFile<BMP>(oxy.createBMP(absSchemeName, false), absFileName);
			writingFile<BMP>(oxy.createBMP(argSchemeName, true), argFileName);
		}
		auto oxz = superposition;
		oxz.transform(a, b, n, u, wavelength, transformType::fractionalFourier, z_begin, z_end, z_n, f);
		absFileName = "oxz_abs.bmp";
		writingFile<BMP>(oxz.createBMP(absSchemeName, false), absFileName);
	}
	catch (runtime_error & e) {
		cerr << endl << "Ошибка! " << e.what() << endl;
		system("pause");
	}
	catch (...) {
		cerr << "Неизвестная ошибка!" << endl;
	}
	return 0;
}