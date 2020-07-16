#include "scheme.h"
#include "CollinsIntegralGPU.h"

scheme::scheme() {
	colorScheme = vector<vector<unsigned char>>();
	schemeName = "";
}

scheme::scheme(string schemeName) {
	vector<pixel> colors;
	colors = loadingData<pixel>(schemeName + ".txt");
	for_each(colors.rbegin(), colors.rend(), [&](pixel pixel) {
		colorScheme.push_back(pixel);
	});
}

scheme::operator vector<vector<unsigned char>>() {
	return colorScheme;
}

pixel::pixel() {
	colors = vector<unsigned char>();
}

pixel::pixel(unsigned char blue, unsigned char green, unsigned char red, unsigned char alpha) {
	colors = { blue, green, red, alpha };
}

pixel::operator vector<unsigned char>() {
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