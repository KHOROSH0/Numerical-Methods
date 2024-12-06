
#include "MPI_ZEIDEL_LU_QR.h"
#include "Matrix_Vector.h"
std::string toString(Vector vector) {
	std::string a;
	for (int i = 0; i < vector.return_size(); i++) {
		std::string p(16, '\0');
		snprintf(&p[0], p.capacity(), "%.10g", vector.return_one_value(i));
		while (!p.empty() && !p.back()) {
			p.pop_back();
		}
		a += p + "\n";
	}
	return a;
}
Vector return_col_with_size(int col, Matrix& A, int size) {
	Vector result{ A.return_n() - size };
	for (int i = 0; i < A.return_n() - size; i++) {
		result.recording_value(i, A.return_value_of_matrix(i + size, col));

	}
	return Vector{ result };
}
Vector return_Basis_with_size(int size) {
	Vector result{ size, };
	for (int i = 0; i < size; i++) {
		if (i == 0) {
			result.recording_value(i, 1);
		}
		else {
			result.recording_value(i, 0);
		}
	}
	return Vector{ result };
}
Matrix multiple_for_matrix(Vector& a, Vector& b) {
	Matrix result{ a.return_size(), b.return_size() };
	for (int i = 0; i < a.return_size(); i++) {
		for (int j = 0; j < a.return_size(); j++) {
			result.recording_value(i, j, a.return_one_value(i) * b.return_one_value(j));
		}
	}
	return Matrix{ result };
}
Matrix create_matrix_with_less_size(Matrix& R, int size) {
	Matrix result{ R.return_m() - size,R.return_m() - size };
	for (int i = 0; i < R.return_m() - size; i++) {
		for (int j = 0; j < R.return_m() - size; j++) {
			result.recording_value(i, j, R.return_value_of_matrix(i + size, j + size));
		}

	}
	return Matrix{ result };
}
Matrix special_function_for_Q(Matrix& Q, int iter) {
	Matrix result{ Q.return_m() + iter, Q.return_m() + iter };
	for (int i = 0; i < Q.return_m() + iter; i++) {
		for (int j = 0; j < Q.return_m() + iter; j++) {
			if (i < iter || j < iter) {
				if (i == j) {
					result.recording_value(i, j, 1);
				}
				else {
					result.recording_value(i, j, 0);
				}
			}
			else {
				result.recording_value(i, j, Q.return_value_of_matrix(i - iter, j - iter));
			}
		}
	}
	return Matrix{ result };
}
Matrix special_fuction_for_R(Matrix& R, Matrix& R_, int iter) {
	Matrix result{ R.return_m(), R.return_m() };
	for (int i = 0; i < R.return_m(); i++) {
		for (int j = 0; j < R.return_m(); j++) {
			if (i < iter || j < iter) {
				result.recording_value(i, j, R.return_value_of_matrix(i, j));
			}
			else {
				result.recording_value(i, j, R_.return_value_of_matrix(i - iter, j - iter));
			}
		}
	}
	return Matrix{ result };
}
int iter_of_Zeidel(Matrix& AA, Vector& bb, double eps) {
	Matrix A{ AA };
	Vector b{ bb };
	Matrix AT{ A.transposed_matrix() };
	double sum = 0;
	for (int i = 0; i < A.return_m(); i++) {
		for (int j = 0; j < A.return_m(); j++) {
			if (i != j) {
				sum = sum + abs(A.return_value_of_matrix(i, j));
			}
		}
	}
	for (int i = 0; i < A.return_m(); i++) {
		if (abs(A.return_value_of_matrix(i, i)) < sum) {
			A = AT * A;
			b = b * AT;
			break;
		}
	}
	sum = 0;
	Vector xkk{ b };
	Matrix C{ A };
	Vector d{ b.return_size() };
	for (int i = 0; i < C.return_m(); i++) {
		d.recording_value(i, b.return_one_value(i) / A.return_value_of_matrix(i, i));
		for (int j = 0; j < C.return_m(); j++) {
			if (i != j) {
				C.recording_value(i, j, ((-1) * A.return_value_of_matrix(i, j) / A.return_value_of_matrix(i, i)));
			}
			else {
				C.recording_value(i, j, 0);
			}
		}
	}
	Vector xk{ d };
	//C.print_matrix();
	//d.print();

	int k = 0;
	while ((xk * A - b).Norma_2() > eps) {
		double step = 0;
		k++;
		for (int i(0); i < A.return_m(); i++) {
			for (int j(0); j < A.return_n(); j++) {
				step = C.return_value_of_matrix(i, j) * xk.return_one_value(j) + step;

			}
			//std::cout << step << std::endl;
			//xk.print();
			xkk.recording_value(i, step + d.return_one_value(i));
			xk = xkk;
			//xk.print();
			step = 0;

		}
	}
	return k;

}
int iter_of_MPI(Matrix AA, Vector& bb, double eps) {
	Matrix A{ AA };
	Vector b{ bb };
	Matrix B{ A };
	Matrix E{ A };
	E.Ed();
	Vector b1{ b };
	Vector c{ b };
	c = ((b1) * (1 / A.Norma_8())); // c = u*b
	//c.print();
	B = (A * (-1 / A.Norma_8())) + E; // B = E - u*A 
	Vector xkk{ c };
	Vector xk{ b };
	Vector dx{ b };
	dx = xkk - xk;
	//dx.print();
	double B_Norma_1 = B.Norma_1();
	double B_Norma_8 = B.Norma_8();
	int k = 0;
	if (B.Norma_1() < 1) {
		while ((B.Norma_1() / (1 - B.Norma_1())) * (dx.Norma_1())) {
			xk = xkk;
			xkk = xk * B + c;
			dx = xkk - xk;
			k++;
		}
	}
	else if (B_Norma_8 < 1) {
		while ((B_Norma_8 / (1 - B_Norma_8)) * (dx.Norma_8()) > eps) {
			xk = xkk;
			xkk = xk * B + c;
			dx = xkk - xk;
			k++;
		}
	}
	else {
		Matrix AT{ A };
		AT = A.transposed_matrix();
		Matrix ATA{ A };
		ATA = AT * A;
		//A.print_matrix();
		Vector b_{ b };
		b_ = b * AT;
		Vector c_{ c };
		c_ = (b_ * (1 / ATA.Norma_8()));
		Matrix B_{ B };
		Vector xkk_{ c_ };
		Vector xk_{ b_ };
		Vector temp{ xkk_ };
		B_ = ATA * (-1 / ATA.Norma_8()) + E;
		temp = xkk_ * ATA - b_;
		//c_.print();
		//b_.print();
		//B_.print_matrix();
		//ATA.print_matrix();
		while (temp.Norma_8() > eps) {
			xk_ = xkk_;
			xkk_ = xk_ * B_ + c_;
			temp = xkk_ * ATA - b_;
			k++;
		}
		return k;
	}
	return k;

}
Vector MPI(Matrix& AA, Vector& bb, double eps) {
	Matrix A{ AA };
	Vector b{ bb };
	Matrix B{ A };
	Matrix E{ A };
	E.Ed();
	Vector b1{ b };
	Vector c{ b };
	c = ((b1) * (1 / A.Norma_8())); // c = u*b
	//c.print();
	B = (A * (-1 / A.Norma_8())) + E; // B = E - u*A 
	Vector xkk{ c };
	Vector xk{ b };
	Vector dx{ b };
	dx = xkk - xk;
	//dx.print();
	double B_Norma_1 = B.Norma_1();
	double B_Norma_8 = B.Norma_8();
	int k = 0;
	if (B.Norma_1() < 1) {
		while ((B.Norma_1() / (1 - B.Norma_1())) * (dx.Norma_1())) {
			xk = xkk;
			xkk = xk * B + c;
			dx = xkk - xk;
			k++;
		}
	}
	else if (B_Norma_8 < 1) {
		while ((B_Norma_8 / (1 - B_Norma_8)) * (dx.Norma_8()) > eps) {
			xk = xkk;
			xkk = xk * B + c;
			dx = xkk - xk;
			k++;
		}
	}
	else {
		Matrix AT{ A };
		AT = A.transposed_matrix();
		Matrix ATA{ A };
		ATA = AT * A;
		//A.print_matrix();
		Vector b_{ b };
		b_ = b * AT;
		Vector c_{ c };
		c_ = (b_ * (1 / ATA.Norma_8()));
		Matrix B_{ B };
		Vector xkk_{ c_ };
		Vector xk_{ b_ };
		Vector temp{ xkk_ };
		B_ = ATA * (-1 / ATA.Norma_8()) + E;
		temp = xkk_ * ATA - b_;
		//c_.print();
		//b_.print();
		//B_.print_matrix();
		//ATA.print_matrix();
		while (temp.Norma_8() > eps) {
			xk_ = xkk_;
			xkk_ = xk_ * B_ + c_;
			temp = xkk_ * ATA - b_;
			k++;
		}
		return xkk_;
	}
	return xkk;

}
Vector Zeidel(Matrix& AA, Vector& bb, double eps) {
	Matrix A{ AA };
	Vector b{ bb };
	Matrix AT{ A.transposed_matrix() };
	double sum = 0;
	for (int i = 0; i < A.return_m(); i++) {
		for (int j = 0; j < A.return_m(); j++) {
			if (i != j) {
				sum = sum + abs(A.return_value_of_matrix(i, j));
			}
		}
	}
	for (int i = 0; i < A.return_m(); i++) {
		if (abs(A.return_value_of_matrix(i, i)) < sum) {
			A = AT * A;
			b = b * AT;
			break;
		}
	}
	sum = 0;
	Vector xkk{ b };
	Matrix C{ A };
	Vector d{ b.return_size() };
	for (int i = 0; i < C.return_m(); i++) {
		d.recording_value(i, b.return_one_value(i) / A.return_value_of_matrix(i, i));
		for (int j = 0; j < C.return_m(); j++) {
			if (i != j) {
				C.recording_value(i, j, ((-1) * A.return_value_of_matrix(i, j) / A.return_value_of_matrix(i, i)));
			}
			else {
				C.recording_value(i, j, 0);
			}
		}
	}
	Vector xk{ d };
	//C.print_matrix();
	//d.print();

	int k = 0;
	while ((xk * A - b).Norma_2() > eps) {
		double step = 0;
		k++;
		for (int i(0); i < A.return_m(); i++) {
			for (int j(0); j < A.return_n(); j++) {
				step = C.return_value_of_matrix(i, j) * xk.return_one_value(j) + step;

			}
			//std::cout << step << std::endl;
			//xk.print();
			xkk.recording_value(i, step + d.return_one_value(i));
			xk.recording_value(i, xkk.return_one_value(i));
			//xk.print();
			step = 0;

		}
	}
	return xk;

}
Vector LU(Matrix& AA, Vector& bb) {
	Matrix A{ AA };
	Vector b{ bb };
	Matrix M{ A };
	Matrix P{ A };
	Matrix E{ A };
	Vector y{ b };
	Vector x{ b };
	P.Ed();
	E.Ed();
	Matrix L{ E };
	Matrix U{ E };
	Matrix Temp{ P };
	int max = 0;
	//M.print_matrix();
	for (int i = 0; i < M.return_m(); i++) {
		for (int j = i; j < M.return_m(); j++) {
			if (abs(M.return_value_of_matrix(j, i)) > max) {
				max = j;
			}
		}
		if (max != 0) {
			P.swap_str(max, i);
			E.swap_str(max, i);
			M = E * M;
			//M.print_matrix();
			E.Ed();
			//M.print_matrix();
			max = 0;
		}
		for (int k(i); k < M.return_m(); k++) {
			for (int t(i); t < M.return_n(); t++) {
				if (k <= i) {
					continue;
				}
				else if (t == i) {
					M.recording_value(k, t, M.return_value_of_matrix(k, t) / M.return_value_of_matrix(i, i));
				}
				else {
					M.recording_value(k, t, M.return_value_of_matrix(k, t) - M.return_value_of_matrix(k, i) * M.return_value_of_matrix(i, t));

				}
			}
			//M.print_matrix();
		}
	}
	for (int i(0); i < M.return_m(); i++) {
		for (int j(0); j < M.return_m(); j++) {
			if (j < i) {
				L.recording_value(i, j, M.return_value_of_matrix(i, j));
			}
			else {
				if (i == j) {
					U.recording_value(i, j, M.return_value_of_matrix(i, j));
				}
				else {
					U.recording_value(i, j, M.return_value_of_matrix(i, j));
				}
			}
		}
	}
	//L.print_matrix();
	//U.print_matrix();
	b = b * P;
	y = b;
	for (int i = 0; i < M.return_m(); i++) { // ищем y
		double step = 0;
		for (int k = 0; k < i; k++) {
			step = step + y.return_one_value(k) * L.return_value_of_matrix(i, k);

		}
		y.recording_value(i, (b.return_one_value(i) - step) / L.return_value_of_matrix(i, i));
	}
	x = y;
	for (int i = M.return_m() - 1; i >= 0; i--) { // ищем x
		double step = 0;
		for (int j = i + 1; j < M.return_m(); j++) {
			step = step + x.return_one_value(j) * U.return_value_of_matrix(i, j);
		}
		x.recording_value(i, (y.return_one_value(i) - step) / U.return_value_of_matrix(i, i));
	}
	return x;


}
Vector QR(Matrix& AA, Vector& bb) {

	Matrix A{ AA };                    // матрица A
	Matrix E{ A };                     // единичная
	E.Ed();
	Vector b{ bb };                    // вектор b
	std::vector <Matrix> Q{ E };         // вектор из ортогоналных матриц Q с начальной E
	std::vector <Matrix> Q_{};         // вектор из ортогональныз матриц меньшей размерности
	std::vector <Matrix> R{ A };       // вектор из матриц которые в конечном счете должны дать матрицу R верхнетреугольную, с начальной A
	std::vector <Matrix> R__{};
	std::vector <Matrix> R_{};         // вектор из R меньшей размерности для получения верхнетреугольной матрицы R 
	std::vector <Vector> w{ };         // вектор из векторов w, тут могут быть вектора разной размерности
	std::vector <Vector> y{};          // вектор из векторов y
	std::vector<Vector> z{};           // базисный вектор 1 0 0...
	std::vector <double> a;            // евклидова норма вектора y
	std::vector <double> r;            // евклидова норма вектора y - z*a



	for (int i = 0; i < A.return_m() - 1; i++) {
		Matrix E{ A.return_m() - i,A.return_m() - i };
		E.Ed();

		y.push_back(return_col_with_size(i, R[i], i));				// создали вектор y из i го столбца матрицы 
		//y[i].print();
		a.push_back(y[i].Norma_2());								// нашли норму y 
		//std::cout << a[i] << std::endl;
		z.push_back(return_Basis_with_size(b.return_size() - i));	//нашли 1 0 0 0... с ь - i - 1 нулем 
		//z[i].print();
		r.push_back((z[i] * a[i] * (-1) + y[i]).Norma_2());			// норма y - a*z
		//std::cout << r[i] << std::endl;
		w.push_back((z[i] * a[i] * (-1) + y[i]) * (1 / r[i]));		//сформировали вектор w
		//w[i].print();
		if (i == 0) {													// так формируем матрицы на первой итерации
			Q.push_back(multiple_for_matrix(w[i], w[i]) * (-2) + E);    //сформировали матрицу Q
			//Q[i + 1].print_matrix();

			R.push_back(Q[i + 1] * R[i]);									// формируес матрицу R
			//R[i + 1].print_matrix(); 
		}
		else if (i >= 1) {															// вторая и последующие итерации (из-за чего у меня все ломается !!!!!!)
			Q_.push_back(multiple_for_matrix(w[i], w[i]) * (-2) + E);
			//Q_[i - 1].print_matrix();
			//R[i].print_matrix();// тут на второй итерации будет 1 матрица в Q_
			R_.push_back(create_matrix_with_less_size(R[i], i)); // теперь у нас на второй итерации будет 2 матрицы в R и одна матрица в R_, следи за этим!!!!!!! То есть разница между их количеством 1!
			//R_[i - 1].print_matrix();
			//(Q_[i - 1] * R_[i - 1]).print_matrix();
			R__.push_back(Q_[i - 1] * R_[i - 1]);
			//R__[i-1].print_matrix();                // теперь в R_ 2 элемента на 2 итерации
			Q.push_back(special_function_for_Q(Q_[i - 1], i)); // делаем матрицу отражений 
			R.push_back(special_fuction_for_R(R[i], R__[i - 1], i));
			//Q[i+1].print_matrix();
			//R[i + 1].print_matrix();
		}
	}
	for (int i = 1; i < A.return_m(); i++) {
		//Q[i].print_matrix();
		Q[i] = Q[i - 1] * Q[i];
	}
	//Q[A.return_m() - 1].print_matrix();
	Vector Y{ b };
	Vector X{ b };
	Matrix QT{ A.return_m(),A.return_m() };
	QT = Q[b.return_size() - 1].transposed_matrix();
	//QT.print_matrix();
	Y = b * QT;
	//Y.print();
	for (int i = A.return_m() - 1; i >= 0; i--) { // ищем x
		double step = 0;
		for (int j = i + 1; j < A.return_m(); j++) {
			step = step + X.return_one_value(j) * R[b.return_size() - 1].return_value_of_matrix(i, j);
		}

		X.recording_value(i, (Y.return_one_value(i) - step) / R[b.return_size() - 1].return_value_of_matrix(i, i));

	}
	return X;
}
