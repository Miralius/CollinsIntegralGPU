#include "Main.h"
#include "ui_main.h"
#include <QtCore/QCoreApplication>

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

//class App : public QMainWindow, Ui_MainWindow
//{
//public:
//	App()
//	{
//		this->setupUi(this);
//	}
//};

int main(int argc, char* argv[]) {

	/*QCoreApplication a(argc, argv);

	auto window = App();
	window.setWindowTitle("ыыы");
	window.show();

	return a.exec();*/
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	try {
		auto a = -2.;
		auto b = 2.;
		auto n = 500;
		auto solitone = field(a, b, n);
		auto sigma = 1.;
		auto eta = 1.;
		auto r = 2.;
		auto number = 3;

		std::vector<double> angleVector = { 0., M_PI / 4., M_PI / 2., M_PI * 3. / 4., M_PI };
		for (auto angle : angleVector) {
			std::string angleName;
			if (angle == 0.) {
				angleName = "0°";
			}
			else if (angle == M_PI / 4.) {
				angleName = "45°";
			}
			else if (angle == M_PI / 2.) {
				angleName = "90°";
			}
			else if (angle == M_PI * 3. / 4.) {
				angleName = "135°";
			}
			else if (angle == M_PI) {
				angleName = "180°";
			}
			std::vector<std::vector<double>> parameterVector = { { 5., 5., -5., -5. }, { 1., 1., 0., 0. }, { 3., 3., 0., 0. }, { 5., 5., 0., 0. }, { 5., 5., 1., 1. }, { 5., 5., 3., 3. }, { 5., 5., 5., 5. } };
			for (auto& parameters : parameterVector) {
				auto alpha = parameters.at(0);
				auto beta = parameters.at(1);
				auto alpha0 = parameters.at(2);
				auto beta0 = parameters.at(3);
				auto parameterName = std::to_string(static_cast<int>(alpha)) + std::to_string(static_cast<int>(beta)) + std::to_string(static_cast<int>(alpha0)) + std::to_string(static_cast<int>(beta0));
				std::cout << "Моделирование входного пучка c параметрами angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
				solitone.airyMode(alpha, beta, alpha0, beta0);
				auto aperture = field(a, b, n);
				aperture.gaussMode(sigma, 1);
				solitone *= aperture;
				solitone.rotate(angle);
				std::string absSchemeName = "fire";
				std::string argSchemeName = "grays";
				auto u = 0.;
				auto wavelength = 650. / 1000000;
				auto z_begin = 100.;
				auto z_end = 1900.;
				auto z_n = 2000;
				auto f = 1000.;
				for (auto i = 0; i <= 2000; i += 50) {
					std::cout << "Моделирование ДрПФ пучка при z = " << i << ", angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
					std::string absFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Solitone\\RotatedBeams\\" + angleName + "\\" + parameterName + "\\" + "oxy_abs_" + std::to_string(i) + ".bmp";
					std::string argFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Solitone\\RotatedBeams\\" + angleName + "\\" + parameterName + "\\" + "oxy_arg_" + std::to_string(i) + ".bmp";
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
				std::cout << "Моделирование продольного сечения пучка c параметрами angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
				field oxz = solitone;
				std::string absFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Solitone\\RotatedBeams\\" + angleName + "\\" + parameterName + "\\" + "oxz_abs.bmp";
				if (std::filesystem::exists(absFileName)) {
					std::cout << "Данный пучок уже смоделирован!" << std::endl;
					continue;
				}
				oxz.ovzFractionalFourierTransform(a, b, n, z_begin, z_end, z_n, wavelength, u, f);
				BMP test = oxz.createBMP(absSchemeName, false);
				writingFile<BMP>(test, absFileName);
			}
		}

		a = -2.;
		b = 2.;
		std::vector<double> ksiVector = { 0., 1., -1., 3., -3. };

		for (auto ksi : ksiVector) {
			auto ksiName = std::to_string(static_cast<int>(ksi));
			std::vector<double> angleVector = { 0., M_PI / 4., M_PI / 2., M_PI * 3. / 4., M_PI };
			for (auto angle : angleVector) {
				std::string angleName;
				if (angle == 0.) {
					angleName = "0°";
				}
				else if (angle == M_PI / 4.) {
					angleName = "45°";
				}
				else if (angle == M_PI / 2.) {
					angleName = "90°";
				}
				else if (angle == M_PI * 3. / 4.) {
					angleName = "135°";
				}
				else if (angle == M_PI) {
					angleName = "180°";
				}
				std::vector<std::vector<double>> parameterVector = { { 5., 5., -5., -5. }, { 1., 1., 0., 0. }, { 3., 3., 0., 0. }, { 5., 5., 0., 0. }, { 5., 5., 1., 1. }, { 5., 5., 3., 3. }, { 5., 5., 5., 5. } };
				for (auto &parameters : parameterVector) {
					auto alpha = parameters.at(0);
					auto beta = parameters.at(1);
					auto alpha0 = parameters.at(2);
					auto beta0 = parameters.at(3);
					auto parameterName = std::to_string(static_cast<int>(alpha)) + std::to_string(static_cast<int>(beta)) + std::to_string(static_cast<int>(alpha0)) + std::to_string(static_cast<int>(beta0));
					std::cout << "Моделирование входного повернутого набора пучков c параметрами ksi = " << static_cast<int>(ksi) << ", angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
					solitone.airyMode(alpha, beta, alpha0, beta0);
					auto aperture = field(a, b, n);
					aperture.gaussMode(sigma, 1);
					solitone *= aperture;
					auto superposition = getRadialSymmetricClusterMode(solitone, sigma, r, number, ksi, eta);
					superposition.rotate(angle);
					std::string absSchemeName = "fire";
					std::string argSchemeName = "grays";
					auto u = 0.;
					auto wavelength = 650. / 1000000;
					auto z_begin = 100.;
					auto z_end = 1900.;
					auto z_n = 2000;
					auto f = 1000.;
					for (auto i = 0; i <= 2000; i += 50) {
						std::cout << "Моделирование ДрПФ повернутого набора пучков при z = " << i << ", ksi = " << static_cast<int>(ksi) << ", angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
						std::string absFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Cluster\\RotatedCluster\\" + ksiName + "\\" + angleName + "\\" + parameterName + "\\" + "oxy_abs_" + std::to_string(i) + ".bmp";
						std::string argFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Cluster\\RotatedCluster\\" + ksiName + "\\" + angleName + "\\" + parameterName + "\\" + "oxy_arg_" + std::to_string(i) + ".bmp";
						if (std::filesystem::exists(absFileName)) {
							std::cout << "Данный повернутый набор пучков уже смоделирован!" << std::endl;
							continue;
						}
						field oxy = superposition;
						oxy.ouvFractionalFourierTransform(a, b, n, wavelength, i, f);
						BMP test = oxy.createBMP(absSchemeName, false);
						BMP test2 = oxy.createBMP(argSchemeName, true);
						writingFile<BMP>(test, absFileName);
						writingFile<BMP>(test2, argFileName);
					}
					std::cout << "Моделирование продольного сечения повернутого набора пучков c параметрами ksi = " << static_cast<int>(ksi) << ", angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
					field oxz = superposition;
					std::string absFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Cluster\\RotatedCluster\\" + ksiName + "\\" + angleName + "\\" + parameterName + "\\" + "oxz_abs.bmp";
					if (std::filesystem::exists(absFileName)) {
						std::cout << "Данный повернутый набор пучков уже смоделирован!" << std::endl;
						continue;
					}
					oxz.ovzFractionalFourierTransform(a, b, n, z_begin, z_end, z_n, wavelength, u, f);
					BMP test = oxz.createBMP(absSchemeName, false);
					writingFile<BMP>(test, absFileName);
				}
			}
		}

		for (auto ksi : ksiVector) {
			auto ksiName = std::to_string(static_cast<int>(ksi));
			std::vector<double> angleVector = { 0., M_PI / 4., M_PI / 2., M_PI * 3. / 4., M_PI };
			for (auto angle : angleVector) {
				std::string angleName;
				if (angle == 0.) {
					angleName = "0°";
				}
				else if (angle == M_PI / 4.) {
					angleName = "45°";
				}
				else if (angle == M_PI / 2.) {
					angleName = "90°";
				}
				else if (angle == M_PI * 3. / 4.) {
					angleName = "135°";
				}
				else if (angle == M_PI) {
					angleName = "180°";
				}
				std::vector<std::vector<double>> parameterVector = { { 5., 5., -5., -5. }, { 1., 1., 0., 0. }, { 3., 3., 0., 0. }, { 5., 5., 0., 0. }, { 5., 5., 1., 1. }, { 5., 5., 3., 3. }, { 5., 5., 5., 5. } };
				for (auto& parameters : parameterVector) {
					auto alpha = parameters.at(0);
					auto beta = parameters.at(1);
					auto alpha0 = parameters.at(2);
					auto beta0 = parameters.at(3);
					auto parameterName = std::to_string(static_cast<int>(alpha)) + std::to_string(static_cast<int>(beta)) + std::to_string(static_cast<int>(alpha0)) + std::to_string(static_cast<int>(beta0));
					std::cout << "Моделирование входного набора повернутых пучков c параметрами ksi = " << static_cast<int>(ksi) << ", angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
					solitone.airyMode(alpha, beta, alpha0, beta0);
					auto aperture = field(a, b, n);
					aperture.gaussMode(sigma, 1);
					solitone *= aperture;
					solitone.rotate(angle);
					auto superposition = getRadialSymmetricClusterMode(solitone, sigma, r, number, ksi, eta);
					std::string absSchemeName = "fire";
					std::string argSchemeName = "grays";
					auto u = 0.;
					auto wavelength = 650. / 1000000;
					auto z_begin = 100.;
					auto z_end = 1900.;
					auto z_n = 2000;
					auto f = 1000.;
					for (auto i = 0; i <= 2000; i += 50) {
						std::cout << "Моделирование ДрПФ набора повернутых пучков при z = " << i << ", ksi = " << static_cast<int>(ksi) << ", angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
						std::string absFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Cluster\\RotatedBeams\\" + ksiName + "\\" + angleName + "\\" + parameterName + "\\" + "oxy_abs_" + std::to_string(i) + ".bmp";
						std::string argFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Cluster\\RotatedBeams\\" + ksiName + "\\" + angleName + "\\" + parameterName + "\\" + "oxy_arg_" + std::to_string(i) + ".bmp";
						if (std::filesystem::exists(absFileName)) {
							std::cout << "Данный набор повернутых пучков уже смоделирован!" << std::endl;
							continue;
						}
						field oxy = superposition;
						oxy.ouvFractionalFourierTransform(a, b, n, wavelength, i, f);
						BMP test = oxy.createBMP(absSchemeName, false);
						BMP test2 = oxy.createBMP(argSchemeName, true);
						writingFile<BMP>(test, absFileName);
						writingFile<BMP>(test2, argFileName);
					}
					std::cout << "Моделирование продольного сечения набора повернутых пучков c параметрами ksi = " << static_cast<int>(ksi) << ", angle = " << angleName << ", alpha|beta|alpha0|beta0 = " << parameterName << std::endl;
					field oxz = superposition;
					std::string absFileName = "C:\\Users\\F-Mir\\OneDrive - ssau.ru\\For hei\\Science work\\2021\\Computer Optics\\Modes\\Cluster\\RotatedBeams\\" + ksiName + "\\" + angleName + "\\" + parameterName + "\\" + "oxz_abs.bmp";
					if (std::filesystem::exists(absFileName)) {
						std::cout << "Данный набор повернутых пучков уже смоделирован!" << std::endl;
						continue;
					}
					oxz.ovzFractionalFourierTransform(a, b, n, z_begin, z_end, z_n, wavelength, u, f);
					BMP test = oxz.createBMP(absSchemeName, false);
					writingFile<BMP>(test, absFileName);
				}
			}
		}

		//auto alpha = 5.;
		//auto beta = 5.;
		//auto alpha0 = -5.;
		//auto beta0 = -5.;
		//auto sigma = 1.;
		//auto ksi = 0.;
		//auto eta = 1.;
		//auto r = 2.;
		//auto number = 3;
		//solitone.airyMode(alpha, beta, alpha0, beta0);
		////solitone.rotate(M_PI);
		//auto aperture = field(a, b, n);
		//aperture.gaussMode(sigma, 1);
		//solitone *= aperture;
		//auto superposition = getRadialSymmetricClusterMode(solitone, sigma, r, number, ksi, eta);
		//superposition.rotate(M_PI);
		//std::string absFileName = "input_abs.bmp";
		//std::string argFileName = "input_phase.bmp";
		//std::string absSchemeName = "fire";
		//std::string argSchemeName = "grays";
		//BMP test = superposition/*solitone*/.createBMP(absSchemeName, false);
		//BMP test2 = superposition/*solitone*/.createBMP(argSchemeName, true);
		//writingFile<BMP>(test, absFileName);
		//writingFile<BMP>(test2, argFileName);
		//auto u = 0.;
		//auto wavelength = 650. / 1000000;
		//auto z_begin = 100.;
		//auto z_end = 1900.;
		//auto z_n = 2000;
		//auto f = 1000.;
		//for (auto i = 250; i <= 1750; i += 250) {
		//	std::cout << std::endl << "Моделирование ДрПФ при z = " << i << std::endl;
		//	field oxy = superposition/*solitone*/;
		//	oxy.ouvFractionalFourierTransform(a, b, n, wavelength, i, f);
		//	absFileName = "oxy_abs_" + std::to_string(i) + ".bmp";
		//	argFileName = "oxy_arg_" + std::to_string(i) + ".bmp";
		//	test = oxy.createBMP(absSchemeName, false);
		//	test2 = oxy.createBMP(argSchemeName, true);
		//	writingFile<BMP>(test, absFileName);
		//	writingFile<BMP>(test2, argFileName);
		//}
		//std::cout << std::endl << "Моделирование продольного сечения пучка" << std::endl;
		//field oxz = superposition/*solitone*/;
		//oxz.ovzFractionalFourierTransform(a, b, n, z_begin, z_end, z_n, wavelength, u, f);
		////oxz.normalize();
		//absFileName = "oxz_abs.bmp";
		//test = oxz.createBMP(absSchemeName, false);
		//writingFile<BMP>(test, absFileName);
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