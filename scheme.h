#ifndef SCHEME_H
#define SCHEME_H

#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
typedef unsigned char byte;

class scheme {
private:
	vector<vector<byte>> colorScheme;
	string schemeName;

public:
	scheme();
	scheme(string schemeName);
	operator vector<vector<byte>>();
};

class pixel {
private:
	vector<byte> colors;

public:
	pixel();
	pixel(byte blue, byte green, byte red, byte alpha);
	operator vector<byte>();
};

istream& operator>>(istream& input, pixel& obj);

#endif