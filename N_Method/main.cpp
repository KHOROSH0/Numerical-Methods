#include "Functions.h"
#include "Matrix_Vector.h"
#include "MPI_ZEIDEL_LU_QR.h"
#include <iostream>
Vector x0{ 2 };
Matrix Y{ 2,2 };
int main() {
	x0* Y;
	double t{ 1e-15 };														
	double eps_for_loc{ 1 };
	double eps_for_N{ 1e-3 };
	Vector ab{ sequential_search(-1 + t, 10, eps_for_loc)};
	double loc1 = ab.return_one_value(0);
	double loc2 = ab.return_one_value(1);
	std::cout << std::endl;
	std::cout << "localisation " << std::endl;
	printf("%10f", loc1);
	printf("%10f", loc2);
	std::cout << std::endl;
	std::cout << "Newton Method" << std::endl;
	printf("%10f", Newton(loc1, loc2, eps_for_N));
	std::cout << std::endl;
	std::cout << "Newton Method 2" << std::endl;
	Newton_for_system(1.3, 0.8, eps_for_N).print();
}