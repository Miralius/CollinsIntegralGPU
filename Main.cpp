#include "Main.h"

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		auto a = 2.;
		auto b = 2.;
		auto n = 100;
		auto output = field(a, b, n);
		auto alpha = 1.;
		auto beta = 1.;
		auto alpha0 = 1.;
		auto beta0 = 1.;
		auto sigma = 1.;
		output.airyMode(alpha, beta, alpha0, beta0, sigma);
		auto aperture = field(a, b, n);
		aperture.gaussMode(sigma, 0);
		output *= aperture;
		string absFileName = "abs.bmp";
		string argFileName = "phase.bmp";
		string absSchemeName = "fire";
		string argSchemeName = "grays";
		
		writingFile<BMP>(output.createBMP(absSchemeName, false), absFileName);
		writingFile<BMP>(output.createBMP(argSchemeName, true), argFileName);
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