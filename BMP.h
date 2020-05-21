#pragma once
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

	static vector<unsigned char> toBinary(vector<int> number) {
		vector<unsigned char> binary;
		for (int i = 0; i < number.at(1); i++) {
			binary.push_back(number.at(0) >> (8 * i));
		}
		return binary;
	}

public:
	BMP();
	BMP(vector<vector<vector<unsigned char>>> picture);
	BMP(const BMP& obj);

	BMP& operator=(const BMP& obj);

	static int toNumber(vector<unsigned char> binary) {
		int temp = 0;
		for (int i = 0; i < binary.size(); i++) {
			temp |= binary.at(i) << (8 * i);
		}
		return temp;
	}

	vector<unsigned char> serialize();
};

istream& operator>>(istream& input, BMP& bmp);
ostream& operator<<(ostream& output, BMP& bmp);

#endif