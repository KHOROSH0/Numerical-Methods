#pragma once
#include "Matrix_Vector.h"
double my_function(double x);												// ������� �� ���������� ����� 1
double derivative(double x);												// ������ �����������									
Vector sequential_search(double a, double b, double eps);					// ���������������� �������
double Newton(double a, double b, double eps);								// ����� �������
Vector Newton_for_system(double a, double b, double eps);

