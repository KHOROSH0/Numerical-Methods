#include "Functions.h"
#include "Matrix_Vector.h"
#include "MPI_ZEIDEL_LU_QR.h"
#include <cmath>

double my_function(double x)
{
	double y = x * log(x + 1) - 0.3;
	return y;
}

double  derivative(double x) {
	double y = (log(x + 1) + x / (x + 1));
	return y;
}
Vector F(double x, double y) {
	Vector result{ 2 };
	result.recording_value(0, sin(x - 1) + y - 1.3);
	result.recording_value(1, x - sin(y+1) - 0.8);
	return result;

}
Matrix Jacobi(double x, double y, double in = 1) {
	Matrix J{ 2,2 };
	J.recording_value(0, 0, cos(x - 1) * in);
	J.recording_value(0, 1, 1);
	J.recording_value(1, 0, 1);
	J.recording_value(1, 1, -cos(y+1) * in);
	return Matrix{ J };
}
Vector sequential_search(double a, double b, double eps)
{
	double N = b - a;
	double k = 1;
	double xkk = b;
	double xk = a;
	double h = (xkk - xk) / N;
	while (abs(xkk - xk) >  eps) {					// критерий остановки
		bool flag = true;
		a = xk;
		b = xkk;
		for (int i = 0; i < N; i++) {
			xk = xk + k * h;
			xkk = xk + k * h;
			
			if (my_function(xk) * my_function(xkk) < 0) {		// если попали в интервал уменьшаем шаг
				h = h / 10;
				flag = false;
				break;
			}
		}
		if (flag == true) {						// если прошли весь интервал и поняли, что ни разу не нашли корень (значит шаг слишком большой)
			h = h / 10;
			xk = a;
			xkk = b;
			flag = false;
		}
	}
	Vector result{ 2 };
	result.recording_value(0, xk);
	result.recording_value(1, xkk);
	return result;
}
double Newton(double a, double b, double eps) {
	double xkk = a;
	double xk = b;
	double ak = a;
	double bk = b;
	double akk = a;
	double bkk = b;
	while (abs(xkk - xk) > 2 * eps) {
		xk = xkk;
		xkk = xk - my_function(xk) / derivative(xk);
		double c = my_function(xkk);
		if (xkk < ak || xkk > bk) {						// как раз если перепрыгнули интервал
			xkk = (ak + bk) / 2;
		}
		if (c > 0) {
			ak = akk;
			akk = ak;
			bkk = xkk;
		}
		if (c <= 0) {
			bk = bkk;
			akk = xkk;
			bkk = bk;
		}
	}

	return xkk;
}
Vector Newton_for_system(double a, double b, double eps) {
	int N = 5;
	Vector x0k{ 2 };
	x0k.recording_value(0, 1.3);					// тут мы заполняем начальное приблежение для начального приближение(круто да?)
	x0k.recording_value(1, 0.8);
	Vector x0kk{ 2 };
	for (int i = 0; i < N; i++) {
		x0k = x0kk;
		Matrix J{ Jacobi(x0k.return_one_value(0), x0k.return_one_value(1), (i) / N) };
		Vector F_{ F(0, 0) };
		x0kk = LU(J, F_);
	}
	Vector xk{ 2 };
	Vector xkk{ x0kk };
	Vector gh{ 2 };
	while (((xkk - xk)).Norma_2() > eps) {
		xk = xkk;
		Matrix J{ Jacobi(xk.return_one_value(0),xk.return_one_value(1)) };
		Vector F_{ F(xk.return_one_value(0), xk.return_one_value(1)) };
		F_ = F_ * (-1);
		Vector gh{ QR(J, F_)};
		xkk = xk + gh;
	}
	return xkk;

}
