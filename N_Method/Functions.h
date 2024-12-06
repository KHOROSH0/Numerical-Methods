#pragma once
#include "Matrix_Vector.h"
double my_function(double x);												// функция из приложения номер 1
double derivative(double x);												// первая производная									
Vector sequential_search(double a, double b, double eps);					// последовательный перебор
double Newton(double a, double b, double eps);								// метод Ньютона
Vector Newton_for_system(double a, double b, double eps);

