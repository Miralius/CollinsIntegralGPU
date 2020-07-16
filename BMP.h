#ifndef BMP_H
#define BMP_H

#include <vector>
#include <Windows.h>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

class BMP {
private:
	static int const countRGBChannel = 4;
	static int const BMPFILEHEADERsize = 14;
	static int const BMPINFOHEADERsize = 124;
	static int const COLORPROFILEsize = 12;
	vector<vector<vector<unsigned char>>> pixels;
	vector<vector<int>> bmpFileHeader; //This vector replaces the struct BITMAPFILEHEADER. Each first element — value, each second element — size of type of var.
	vector<vector<int>> bmpInfoHeader; //This vector replaces the struct BITMAPINFOHEADER. Each first element — value, each second element — size of type of var.
	vector<vector<int>> colorProfile;

	void initHeaders(int width, int height);

public:
	BMP();
	BMP(vector<vector<vector<unsigned char>>> picture);
	BMP(const BMP& obj);
	static vector<unsigned char> toBinary(vector<int> number);
	static int toNumber(vector<unsigned char> binary);
	BMP& operator=(const BMP& obj);
	operator vector<unsigned char>();
};

istream& operator>>(istream& input, BMP& bmp);
ostream& operator<<(ostream& output, BMP& bmp);

#endif