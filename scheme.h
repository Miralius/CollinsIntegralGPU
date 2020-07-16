#ifndef SCHEME_H
#define SCHEME_H

#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

class scheme {
private:
	vector<vector<unsigned char>> colorScheme;
	string schemeName;

public:
	scheme();
	scheme(string schemeName);
	operator vector<vector<unsigned char>>();
};

class pixel {
private:
	vector<unsigned char> colors;

public:
	pixel();
	pixel(unsigned char blue, unsigned char green, unsigned char red, unsigned char alpha);
	operator vector<unsigned char>();
};

istream& operator>>(istream& input, pixel& obj);

#endif