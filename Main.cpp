#include "Main.h"

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		auto a = 2.;
		auto b = 2.;
		auto n = 500;
		auto output = field(a, b, n);
		auto alpha = 1.;
		auto beta = 1.;
		auto alpha0 = 1.;
		auto beta0 = 1.;
		auto sigma = 1.;
		output.airyMode(alpha, beta, alpha0, beta0, sigma);
		auto aperture = field(a, b, n);
		aperture.gaussMode(sigma, 0.);
		output *= aperture;
		auto u = 0.;
		auto wavelength = 650. / 1000000;
		auto z = 1000.;
		auto z_begin = 0.001;
		auto z_end = 1000.;
		auto z_n = 1000;
		auto f = 1000.;
		auto oxy = output;
		oxy.transform(a, b, n, wavelength, transformType::fractionalFourier, z, f);
		auto oxz = output;
		oxz.transform(a, b, n, u, wavelength, transformType::fractionalFourier, z_begin, z_end, z_n, f);
		oxz.normalize();
		string absFileName = "input_abs.bmp";
		string argFileName = "input_phase.bmp";
		string absSchemeName = "fire";
		string argSchemeName = "grays";
		writingFile<BMP>(output.createBMP(absSchemeName, false), absFileName);
		writingFile<BMP>(output.createBMP(argSchemeName, true), argFileName);

		absFileName = "oxy_abs.bmp";
		writingFile<BMP>(oxy.createBMP(absSchemeName, false), absFileName);
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