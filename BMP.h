#ifndef BMP_H
#define BMP_H

#define _HAS_STD_BYTE 0
#include <vector>
#include <Windows.h>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

class BMP {
private:
	static auto const countRGBChannel = 4;
	static auto const BMPFILEHEADERsize = 14;
	static auto const BMPINFOHEADERsize = 124;
	static auto const COLORPROFILEsize = 12;
	vector<vector<vector<byte>>> pixels;
	vector<vector<int>> bmpFileHeader; //This vector replaces the struct BITMAPFILEHEADER. Each first element — value, each second element — size of type of var.
	vector<vector<int>> bmpInfoHeader; //This vector replaces the struct BITMAPINFOHEADER. Each first element — value, each second element — size of type of var.
	vector<vector<int>> colorProfile;
	void initHeaders(int width, int height);

public:
	BMP();
	BMP(vector<vector<vector<byte>>>& picture);
	BMP(const BMP& obj);
	static vector<byte> toBinary(vector<int> number);
	static int toNumber(vector<byte> binary);
	BMP& operator=(const BMP& obj);
	operator vector<byte>();
};

istream& operator>>(istream& input, BMP& bmp);
ostream& operator<<(ostream& output, BMP& bmp);

#endif