
#include "Matrix_Vector.h"
#include <iostream>
#include <istream>
#include <vector>
#include <string>
Matrix::Matrix(int n, int m, std::istream& sth) {// create
	if (n == 0 || m == 0) {
		get_size(sth);
		get_matrix(sth);
	}
	else if (n != 0 && m != 0) {
		this->n = n;
		this->m = m;
		create_memory(n, m);

	}
}
Matrix::~Matrix() { // delete
	delete_memory();
}
Matrix::Matrix(const Matrix& matrix) { // copy
	n = matrix.n;
	m = matrix.n;
	create_memory(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			mat[i][j] = matrix.mat[i][j];
		}
	}
}
Matrix Matrix::Ed() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				mat[i][j] = 1;
			}
			else {
				mat[i][j] = 0;
			}
		}
	}
	return *this;
}
void Matrix::get_size(std::istream& sth) {
	int t, k;
	sth >> t >> k;
	n = t;
	m = k;
	//delete_memory();
	create_memory(n, m);

}
double Matrix::Norma_1() {
	double sum = 0;
	double max = 0;
	for (int i(0); i < n; i++) {
		for (int j(0); j < m; j++) {
			sum = sum + (abs)(mat[i][j]);
		}
		if (sum > max) {
			max = sum;
		}
		sum = 0;
	}
	return max;
}
double Matrix::Norma_8() {
	double sum = 0;
	double max = 0;
	for (int i(0); i < m; i++) {
		for (int j(0); j < n; j++) {
			sum = sum + (abs)(mat[j][i]);
		}
		if (sum > max) {
			max = sum;
		}
		sum = 0;
	}
	return max;
}
Matrix Matrix:: operator = (const Matrix matrix) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			mat[i][j] = matrix.mat[i][j];
		}
	}
	return   *this;
}
Matrix Matrix:: operator + (Matrix& matr) { // sum
	Matrix result{ n,m };
	if (matr.n == n && matr.m == m) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				result.mat[i][j] = matr.mat[i][j] + mat[i][j];
			}
		}
	}
	return  result;
}
Matrix Matrix:: operator - (Matrix& matrix) { // dif
	Matrix result{ *this };
	if (matrix.n == n && matrix.m == m) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				result.mat[i][j] = mat[i][j] - matrix.mat[i][j];
			}
		}
	}
	return result;
}
Matrix Matrix:: operator * (Matrix& matrix) {
	Matrix result{ n, matrix.m };// multiply matrix and matrix
	if (m == matrix.n) {
		double step = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < matrix.m; j++) {
				for (int s = 0; s < m; s++) {
					step = mat[i][s] * matrix.mat[s][j] + step;
				}
				result.mat[i][j] = step;
				step = 0;
			}
		}
	}
	return  Matrix{ result };
}
Matrix Matrix:: operator * (const double number) {
	// multiply matrix and number
	Matrix result{ n,m };
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			result.mat[i][j] = number * mat[i][j];
		}
	}
	return  result;
}

bool Matrix:: operator == (const Matrix& matrix) { // compare 
	if (matrix.n == n && matrix.m == m) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (mat[i][j] != matrix.mat[i][j]) {
					std::cout << " NO " << std::endl;
					return false;
				}
			}
		}
		std::cout << "YES" << std::endl;
		return true;
	}
	else {
		std::cout << "comparison is impossible" << std::endl;
		return false;
	}
}
Matrix Matrix::transposed_matrix() {
	// here we tronspose the matrix
	Matrix result{ *this };
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < m; j++) {
			double temp = result.mat[i][j];
			result.mat[i][j] = result.mat[j][i];
			result.mat[j][i] = temp;
		}
	}
	return  result;
}
void Matrix::swap_str(int k, int r) {
	for (int i = 0; i < m; i++) {
		std::swap(mat[r][i], mat[k][i]);
	}
}

int Matrix::return_n() {
	return n;
}
double Matrix::return_value_of_matrix(int str, int col) {
	return mat[str][col];
}
void Matrix::recording_value(int str, int col, double value) {
	mat[str][col] = value;
}
void Matrix::get_matrix(std::istream& sth) {

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			sth >> mat[i][j];
		}
	}
}
void Matrix::print_matrix() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			std::cout << mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
}
void Matrix::create_memory(int n, int m) {
	mat = new double* [n];
	mat[0] = new double[n * m];
	for (int i = 1; i < n; i++) {
		mat[i] = mat[i - 1] + m;
	}
}
void Matrix::delete_memory() {
	delete mat[0];
	delete[] mat;
}
int Matrix::return_m() {
	return m;
}
Vector::Vector(int size, std::istream& sth) {
	if (size == -1) {
		get_size(sth);
		get_vector(sth);
	}
	else if (size > -1) {
		this->size = size;
		create_memory(size);
		for (int i = 0; i < size; i++) {
			vec[i] = 0;
		}
	}
}
double Vector::Norma_1() {
	double sum = 0;
	for (int i(0); i < size; i++) {
		sum = sum + abs(vec[i]);
	}
	return sum;
}
double Vector::Norma_2() {
	double sum = 0;
	for (int i(0); i < size; i++) {
		sum = sum + (vec[i]) * (vec[i]);
	}
	return sqrt(sum);
}
double Vector::Norma_8() {
	int current = 0;
	for (int i(0); i < size; i++) {
		if (abs(vec[i]) > abs(vec[current])) {
			current = i;
		}
	}
	return (abs)(vec[current]);

}

void Vector::print() {
	for (int i = 0; i < size; i++) {
		std::cout << vec[i] << " ";
	}
	std::cout << std::endl;
}
Vector::Vector(const Vector& vector) {
	size = vector.size;
	create_memory(size);
	for (int i = 0; i < size; i++) {
		vec[i] = vector.vec[i];
	}
}
Vector::~Vector() {
	delete_memory();
}
Vector Vector:: operator * (Matrix& matrix) { // multiply vector and matrix
	Vector result{ size };
	double step = 0;
	for (int i = 0; i < matrix.return_n(); i++) {
		for (int j = 0; j < matrix.return_n(); j++) {
			step = vec[j] * matrix.return_value_of_matrix(i, j) + step;
		}
		result.recording_value(i, step);
		step = 0;
	}
	return  result;
}
bool Vector:: operator == (Vector& v) // compare
{
	if (size != v.size) {
		std::cout << " these vectors are not comparable " << std::endl;
		return false;
	}
	else {
		bool flag = true;
		for (int i = 0; i < size; i++) {
			if (vec[i] != v.vec[i]) {
				return vec[i] == v.vec[i];
			}
		}
		return true;
	}
}
Vector Vector:: operator + (Vector& v) { // sum
	Vector result{ size };
	for (int i = 0; i < size; i++) {
		result.vec[i] = v.vec[i] + vec[i];
	}
	return result;
}
Vector Vector:: operator - (Vector& v) { // dif
	Vector result{ size };
	for (int i = 0; i < size; i++) {
		result.vec[i] = vec[i] - v.vec[i];
	}
	return result;
}
double Vector:: operator * (Vector& v) { // multiply vector and vector
	if (size != v.size) {
		std::cout << " product is impossible " << std::endl;
		return -1;
	}
	else {
		double step = 0;
		for (int i = 0; i < size; i++) {
			step = v.vec[i] * vec[i] + step;
		}
		return step;
	}
}
Vector Vector:: operator = (const Vector vector) {
	for (int i = 0; i < size; i++) {
		vec[i] = vector.vec[i];
	}
	return  *this;
}
Vector Vector:: operator * (double number) {
	Vector result{ *this };// multiply vector and number
	for (int i = 0; i < size; i++) {
		result.vec[i] = vec[i] * number;
	}
	return  result;
}

void Vector::get_size(std::istream& sth) {
	sth >> size;
	create_memory(size);
}
void Vector::get_vector(std::istream& sth) {
	for (int i = 0; i < size; i++) {
		sth >> vec[i];
	}
}

void Vector::swap_vector(int k, int r) {
	std::swap(vec[k], vec[r]);
}
void Vector::change_value_vector(int index_current, int index_dif, double mat_dif, double mat_current) {
	vec[index_current] = vec[index_current] * mat_dif - vec[index_dif] * mat_current;
}
void Vector::value_for_solve(int index_current, double mat_current) {
	vec[index_current] = mat_current;
}
void Vector::solve_for_vector(int  index_current, double mat_current) {
	vec[index_current] = vec[index_current] / mat_current;
}
int Vector::return_size() {
	return size;
}
double Vector::return_values() {
	return *vec;
}
double Vector::return_one_value(int i) {
	return vec[i];
}
void Vector::recording_value(int i, double value) {
	vec[i] = value;
}
void Vector::create_memory(int size) {
	vec = new double[size];
}
void Vector::delete_memory() {
	delete[] vec;
	vec = nullptr;
}