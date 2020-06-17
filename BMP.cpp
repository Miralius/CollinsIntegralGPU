#include "BMP.h"

void BMP::initHeaders(int width, int height) {
	bmpFileHeader = { {0x4D42, 2}, {width * height * countRGBChannel + BMPFILEHEADERsize + BMPINFOHEADERsize + COLORPROFILEsize, 4}, {0, 2}, {0, 2}, {BMPFILEHEADERsize + BMPINFOHEADERsize + COLORPROFILEsize, 4} };
	bmpInfoHeader = { {BMPINFOHEADERsize, 4}, {width, 4}, {height, 4}, {1, 2}, {countRGBChannel * 8, 2}, {BI_BITFIELDS, 4}, {0, 4}, {0xEC4, 4}, {0xEC4, 4}, {0, 4}, {0, 4},
		{0x00FF0000, 4}, {0x0000FF00, 4}, {0x000000FF, 4}, {static_cast<int>(0xFF000000)/*-16777216*/, 4}, {LCS_WINDOWS_COLOR_SPACE, 4}, {0, 36}, {0, 4}, {0, 4}, {0, 4},
		{0, 4}, {0, 4}, {0, 4}, {0, 4} };
	colorProfile = { {0x000000FF, 4}, {0x0000FF00, 4}, {0x00FF0000, 4} };
}

BMP::BMP() {
	initHeaders(0, 0);
}

BMP::BMP(vector<vector<vector<unsigned char>>> picture) : pixels(picture) {
	initHeaders((int)picture.size(), (int)picture.at(0).size());
}

BMP::BMP(const BMP& obj) {
	pixels = obj.pixels;
	bmpFileHeader = obj.bmpFileHeader;
	bmpInfoHeader = obj.bmpInfoHeader;
	colorProfile = obj.colorProfile;
}

BMP& BMP::operator=(const BMP& obj) {
	if (this == &obj) return *this;
	pixels = obj.pixels;
	bmpFileHeader = obj.bmpFileHeader;
	bmpInfoHeader = obj.bmpInfoHeader;
	colorProfile = obj.colorProfile;
	return *this;
}

vector<unsigned char> BMP::serialize() {
	vector<unsigned char> serializedBMP;
	serializedBMP.reserve(bmpFileHeader.at(1).at(0));
	for (vector<int> data : bmpFileHeader) {
		for (unsigned char byte : toBinary(data)) {
			serializedBMP.push_back(byte);
		}
	}
	for (vector<int> data : bmpInfoHeader) {
		for (unsigned char byte : toBinary(data)) {
			serializedBMP.push_back(byte);
		}
	}
	for (vector<int> data : colorProfile) {
		for (unsigned char byte : toBinary(data)) {
			serializedBMP.push_back(byte);
		}
	}
	for_each(pixels.rbegin(), pixels.rend(), [&](vector<vector<unsigned char>> row) {
		for (vector<unsigned char> pixel : row) {
			for (unsigned char color : pixel) {
				serializedBMP.push_back(color);
			}
		}
		});
	return serializedBMP;
}

ostream& operator<<(ostream& output, BMP& bmp) {
	vector<unsigned char> data = bmp.serialize();
	for (unsigned char value : data) {
		output << value;
	}
	return output;
}

istream& operator>>(istream& input, BMP& bmp) {
	int const bytesBeforeOffset = 10; //count of bytes which following before pixel's offset
	unsigned char buffer;
	for (int i = 0; i < bytesBeforeOffset; i++) {
		if (!(input >> buffer)) return input;
	}

	vector<unsigned char> bufferVector;
	for (int i = 0; i < sizeof(int); i++) {
		if (!(input >> buffer)) return input;
		bufferVector.push_back(buffer);
	}
	int offset = BMP().toNumber(bufferVector);

	int const bytesAfterOffsetBeforeWidth = 4; //count of bytes between offset and width
	for (int i = 0; i < bytesAfterOffsetBeforeWidth; i++) {
		if (!(input >> buffer)) return input;
	}

	vector<int> size;
	bufferVector.clear();
	for (int i = 0; i < 2; i++) { //2 because width & height
		for (int j = 0; j < sizeof(int); j++) {
			if (!(input >> buffer)) return input;
			bufferVector.push_back(buffer);
		}
		size.push_back(BMP().toNumber(bufferVector));
		bufferVector.clear();
	}

	int bytesAfterSizeBeforeBitCount = 2; //count of bytes between height and bitCount
	for (int i = 0; i < bytesAfterSizeBeforeBitCount; i++) {
		if (!(input >> buffer)) return input;
	}

	for (int i = 0; i < sizeof(short); i++) {
		if (!(input >> buffer)) return input;
		bufferVector.push_back(buffer);
	}
	int bitCount = BMP().toNumber(bufferVector);

	int bytesAfterBitCount = offset - 30; //count of bytes after bitCount
	for (int i = 0; i < bytesAfterBitCount; i++) {
		if (!(input >> buffer)) return input;
	}

	vector<vector<vector<unsigned char>>> pixels;
	vector<unsigned char> pixel;
	for (int i = 0; i < size.at(0); i++) {
		pixels.push_back(vector<vector<unsigned char>>());
		for (int j = 0; j < size.at(1); j++) {
			if (bitCount == 32) {
				for (int k = 0; k < 4; k++) { //because BGR + alpha
					if (!(input >> buffer)) return input;
					pixel.push_back(buffer);
				}
			}
			if (bitCount == 24) {
				for (int k = 0; k < 3; k++) { //because BGR
					if (!(input >> buffer)) return input;
					pixel.push_back(buffer);
				}
				pixel.push_back(255);
			}
			pixels.back().push_back(pixel);
			pixel.clear();
		}
	}
	reverse(pixels.begin(), pixels.end());
	bmp = BMP(pixels);
	return input;
}
