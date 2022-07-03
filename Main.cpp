#include "Main.h"

field getRadialSymmetricClusterMode(const field& solitone, double sigma, double radius, int number, double ksi, double eta)
{
	field superposition = solitone;
	superposition.shift(radius, 0);
	auto parameterMode = field(superposition.getY().back(), -superposition.getY().back(), static_cast<int>(superposition.getY().size()));
	parameterMode.initialTransverseVelocityAndPowerFactorExpMode(ksi, eta, sigma, radius, 0);
	parameterMode.transpose();
	superposition *= parameterMode;
	for (auto i = 1; i < number; i++) {
		field solitoneShifted = solitone;
		auto phi_n = M_PI * 2 * i / number;
		auto c_xn = radius * cos(phi_n);
		auto c_yn = radius * sin(phi_n);
		solitoneShifted.shift(c_xn, c_yn);
		parameterMode.initialTransverseVelocityAndPowerFactorExpMode(ksi, eta, sigma, c_xn, c_yn);
		parameterMode.transpose();
		solitoneShifted *= parameterMode;
		superposition += solitoneShifted;
	}
	return superposition;
}

int main(int argc, char* argv[]) {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		auto a = -4.;
		auto b = 4.;
		auto n = 500;
		auto solitone = field(a, b, n);
		auto sigma = 1.5;
		
		std::cout << "Моделирование входного пучка " << std::endl;
		solitone.doubledLaguerreMode(sigma);
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
			std::string absFileName = "Modes\\Solitone\\oxy_abs_" + std::to_string(i) + ".bmp";
			std::string argFileName = "Modes\\Solitone\\oxy_arg_" + std::to_string(i) + ".bmp";
			if (std::filesystem::exists(absFileName)) {
				std::cout << "Данный пучок уже смоделирован!" << std::endl;
				continue;
			}
			field oxy = solitone;
			oxy.ouvFractionalFourierTransform(a, b, n, wavelength, i, f);
			BMP test = oxy.createBMP(absSchemeName, false);
			BMP test2 = oxy.createBMP(argSchemeName, true);
			writingFile<BMP>(test, absFileName);
			writingFile<BMP>(test2, argFileName);
		}
		std::cout << "Моделирование продольного сечения пучка" << std::endl;
		field oxz = solitone;
		std::string absFileName = "Modes\\Solitone\\oxz_abs.bmp";
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