#include <iostream>
#include "Header.h"
#include <iomanip>
#include "manip.cpp"

double Heron(double x, const double& eps);
int fact(int n);
double sin2(const double& a, double& x, const double& b, const double& eps);
double arctan(double& ResultofHeron, const double& eps);
int main() {
	vector <string> colums = {"X", "U", "Delta_U", "U_Precision", "U_Absolute_mistake", "V", "Delta_V", "V_Precision", "V_Absolute_mistake",
		"F", "Delta_F", "F_Precision", "F_Absolute_mistake","Z", "Delta_Z", "Z_Precision","Z_Absolute_mistake" };
	Writer writer(colums);
	for (int i(0); i < 11; i++) {
		double x1 = 20 + i;
		double x = x1 / 100;
		double sqrt_ = sqrt(0.9 * x + 1)/(1 - x*x);
		double sq = Heron(0.9* x+1, 0.5e-6) / (1 - x * x);
		vector <long double> valume = {x,arctan(sq, 0.5 * 0.5e-6),0.5 * 0.5e-6,  atan(sq), abs((atan(sq) - arctan(sq, 0.5 * 0.5e-6))),
			sin2(3, x, 0.6, 0.5e-6),0,5*1e-6, sin(3 * x + 0.6),
			abs(sin(3 * x + 0.6) - sin2(3, x, 0.6, 0.5e-6)), sq,
			0.5e-6,sqrt_,abs((sqrt_) - sq),
			arctan(sq, 1e-6) - sin2(3, x, 0.6, 0.5e-6), 0.5e-6 , atan(sqrt_) - sin(3 * x + 0.6),
			abs((arctan(sq, 1e-6) - sin2(3, x, 0.6, 0.5e-6)) - (atan(sqrt_) - sin(3*x+0.6)))};
		vector<string> value;

		for (auto val : valume) {
			string s(14, ' ');
			snprintf(&s[0], s.capacity(), "%.14lf", val);
			value.push_back(s);
		}
		writer.add(value);
	
	}
	writer.print();
}