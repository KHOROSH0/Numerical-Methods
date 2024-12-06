#pragma once
#include "Matrix_Vector.h"
Vector MPI(Matrix& AA, Vector& bb, double eps);
Vector Zeidel(Matrix& AA, Vector& bb, double eps);
Vector LU(Matrix& AA, Vector& bb);
Vector QR(Matrix& AA, Vector& bb);
std::string toString(Vector vector);
int iter_of_Zeidel(Matrix& AA, Vector& bb, double eps);
int iter_of_MPI(Matrix AA, Vector& bb, double eps);