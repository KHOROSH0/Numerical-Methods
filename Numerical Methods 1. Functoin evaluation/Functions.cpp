#include <cmath>
#include "Header.h"
const double PI = 3.1415926535;
 double Heron(double x, const double& eps) {
	 double ri = 0;
	 double rii = 2;
	while (abs(rii - ri) > eps) { // sqrt
		ri = rii;
		rii = (double)(ri + ((x) / ri)) / 2;
	}
	return rii;
}

 int fact(int n ) {
	 if (n == 1) {
		 return 1;
	 }
	 else if (n == 0) {
		 return 1;
	 }
		 return (fact(n - 1) * n);
 }
 double sin2(const double& a, double& x, const double& b, const double& eps) {
	 double sumi = -1;
	 double sumii = 0;
	 unsigned int k = 0;
	 if (a * x + b < PI / 4) {
		 while (abs(sumii - sumi) > eps) {
			 sumi = sumii;
			 sumii = sumii + (pow((-1), k) * (pow((a*x + b), 2 * k + 1)) / fact(2 * k + 1));
			 k++;
		 }
	 }
	 else {
		 while (abs(sumii - sumi) > eps) {
			 sumi = sumii;
			 sumii = sumi + (pow((-1), k) * (pow(PI/2 - (a*x+b), 2 * k)) / fact(2 * k));
			 k++;
		 }
	 }
	 return sumii;

 }
 double arctan(double& ResultofHeron,const double& eps) {
	 double sumi = -1;
	 double sumii = 0;
	 unsigned int k = 0;
	 if (abs(ResultofHeron) < 1) {
		 while (abs(sumii - sumi) > eps) {
			 sumi = sumii;
			 sumii = sumii + (pow(-1, k) * pow(ResultofHeron, 2 * k + 1)) / (2 * k + 1);
			 k++;
		 }
		 return sumii;
	 }
	 else if (abs(ResultofHeron) >= 1)
		 while (abs(sumii - sumi) > eps) {
			 sumi = sumii;
			 sumii = (sumii + (pow(-1, k) / (pow(ResultofHeron, 2 * k + 1) * (2 * k + 1))));
			 k++;
		 }
	 return (((PI / 2) * (ResultofHeron / abs(ResultofHeron))) - sumii);
 }