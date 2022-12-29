#pragma once

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include "stdio.h"

#define PI 3.1415926
using namespace std;

class CycloidGear
{
public:
	std::pair<double, double> cycloid_th(double t);
	std::pair<double, double> cycloid(double t);
	std::pair<double, double> cycloid_md(double t);

	CycloidGear() = default;

	CycloidGear(double R0, double r0, double e0, int Z0, double rm0) : R(R0), r(r0), e(e0), Z(Z0), rm(rm0) {}

	// Calculate the Md according to the input bar position on the axis X.
	double printMd(double md0) {
		double md = (md0) * sqrt(2+2*cos(PI/Z)) + 2 * rm;
		printf("Md: %.5f\n",md);

		return md;
	}

	// Generate a point cloud for the curve.
	void curve_output(const std::string& filename, const int curve = 2) {
		std::ofstream ofs(filename);
		std::pair<double, double> (CycloidGear::*f)(double t) = &CycloidGear::cycloid_md;
		switch (curve)
		{
		case 0:
			f = &CycloidGear::cycloid_th;
			break;
		case 1:
			f = &CycloidGear::cycloid;
		default:
			break;
		}
		for (double t = 0; t < 2 * PI; t += 0.001) {
			std::pair<double, double> pol = (this->*f)(t);
			double x = pol.first;
			double y = pol.second;
			ofs << x << "," << y << ",0" << std::endl;
		}
	}

	// Find the position of the bar on the axis X.
	// Return the coordinate (x, y) of bar's center.
	std::pair<double, double> find_bar() {
		// Resolution:
		double dt = 1E-9;
		double threshold = dt * 100;
		// find all point whose y coordinate is zero (or less than 1e-6)
		double t = 0;
		std::pair<double, double> point = cycloid_md(0);
		std::vector<std::pair<double, double>> result;
		double y = point.second;
		while (result.size() <= 2 && t < PI / 78) {
			t += dt;
			point = cycloid_md(t);
			y = point.second;
			if (fabs(y) < threshold)
				result.push_back(point);
		}
		if (result.size() < 2) {
			std::cout << "The position of the bar is not found!" << std::endl;
			return std::pair<double, double>(0, 0);
		}
		else
			return result.back();
	}

	// Return the distance between p1 and p2.
	double distance(std::pair<double, double> p1, std::pair<double, double>p2) {
		return sqrt((p1.first - p2.first) * (p1.first - p2.first) + (p1.second - p2.second) * (p1.second - p2.second));
	}

	void setR(double D) { R = D / 2; }
	void setr(double d) { r = d / 2; }
	void sete(double e) { (*this).e = e; }
	void setZ(int Z) { (*this).Z = Z; }
	void setrm(double dm) { rm = dm / 2; }


private:
	double R;		// Õë³Ý·Ö²¼Ô²°ë¾¶
	double r;		// Õë³Ý°ë¾¶
	double e;		// Æ«ÐÄ¾à
	double rm;		// Á¿°ô°ë¾¶
	int Z;			// °ÚÏßÂÖ³ÝÊý
	vector<std::pair<std::pair<double, double>, double>> vdata; // <<dr, dR>, md>
};

// Calculate the coordinate of point on the cycloid curve according to the input angle.
inline std::pair<double, double> CycloidGear::cycloid_th(double t) {
	double x = 0.0, y = 0.0;
	double Z2 = Z + 1;
	double K = e * Z2 / R;
	double S = 1 + K * K - 2 * K * cos(Z * t);
	double sgamma = (-K * cos(Z2 * t) + cos(t)) / sqrt(S);
	double cgamma = (K * sin(Z2 * t) - sin(t)) / sqrt(S);
	double x0 = R * cos(t) - e * cos(Z2 * t);
	double y0 = R * sin(t) - e * sin(Z2 * t);
	return std::pair<double, double>(x0, y0);
}

// Calculate the coordinate of point on the cycloid curve according to the input angle.
inline std::pair<double, double> CycloidGear::cycloid(double t) {
	double x = 0.0, y = 0.0;
	double Z2 = Z + 1;
	double K = e * Z2 / R;
	double S = 1 + K * K - 2 * K * cos(Z * t);
	double sgamma = (-K * cos(Z2 * t) + cos(t)) / sqrt(S);
	double cgamma = (K * sin(Z2 * t) - sin(t)) / sqrt(S);
	double x0 = R * cos(t) - e * cos(Z2 * t);
	double y0 = R * sin(t) - e * sin(Z2 * t);
	x = x0 - r * sgamma;
	y = y0 + r * cgamma;
	return std::pair<double, double>(x, y);
}

// Calculate the coordinate of bar's center rolling on the cycloid curve according to the input angle.
inline std::pair<double, double> CycloidGear::cycloid_md(double t) {
	double x = 0.0, y = 0.0;
	double Z2 = Z + 1;
	double K = e * Z2 / R;
	double S = 1 + K * K - 2 * K * cos(Z * t);
	double sgamma = (-K * cos(Z2 * t) + cos(t)) / sqrt(S);
	double cgamma = (K * sin(Z2 * t) - sin(t)) / sqrt(S);
	double x0 = R * cos(t) - e * cos(Z2 * t);
	double y0 = R * sin(t) - e * sin(Z2 * t);
	x = x0 - (r - rm) * sgamma; // ¼ÆËãÁ¿°ô¹ö¶¯¹ì¼£
	y = y0 + (r - rm) * cgamma; // ¼ÆËãÁ¿°ô¹ö¶¯¹ì¼£
	return std::pair<double, double>(x, y);
}

