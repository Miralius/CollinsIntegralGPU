#include "scheme.h"
#include "Main.h"

scheme::scheme() {
	colorScheme = vector<vector<byte>>();
	schemeName = "";
}

scheme::scheme(string schemeName) {
	auto colors = loadingData<pixel>("schemes\\" + schemeName + ".txt");
	for_each(colors.rbegin(), colors.rend(), [&](pixel pixel) {
		colorScheme.push_back(pixel);
	});
}

scheme::operator vector<vector<byte>>() {
	return colorScheme;
}

pixel::pixel() : colors(vector<byte>()) {
}

pixel::pixel(byte blue, byte green, byte red, byte alpha) : colors({ blue, green, red, alpha }) {
}

pixel::operator vector<byte>() {
	return colors;
}

istream& operator>>(istream& input, pixel& obj) {
	int blue, green, red;
	input >> red >> green >> blue;
	if (!input) {
		return input;
	}
	obj = pixel(blue, green, red, 255);
	return input;
}