#include "scheme.h"

scheme::scheme() {
	colorScheme = std::vector<std::vector<byte>>();
	schemeName = "";
}

scheme::scheme(std::string schemeName) {
	auto colors = loadingData<pixel>("schemes\\" + schemeName + ".txt");
	colorScheme.reserve(colors.size());
	std::for_each(colors.rbegin(), colors.rend(), [&](pixel pixel) {
		colorScheme.emplace_back(pixel);
	});
}

scheme::operator std::vector<std::vector<byte>>() {
	return colorScheme;
}

pixel::pixel() : colors(std::vector<byte>()) {
}

pixel::pixel(byte blue, byte green, byte red, byte alpha) : colors({ blue, green, red, alpha }) {
}

pixel::operator std::vector<byte>() {
	return colors;
}

std::istream& operator>>(std::istream& input, pixel& obj) {
	int blue, green, red;
	input >> red >> green >> blue;
	if (!input) {
		return input;
	}
	obj = pixel(blue, green, red, 255);
	return input;
}