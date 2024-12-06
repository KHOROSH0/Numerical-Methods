#pragma once
#pragma once
#include <iostream>
#include <istream>
#include <vector>
#include <string>

class Matrix {
	int n, m;
	double** mat;

public:
	Matrix(int n = 0, int m = 0, std::istream& sth = std::cin);
	~Matrix();
	Matrix(const Matrix& matrix);
	void get_size(std::istream& sth);
	double Norma_1();
	double Norma_8();
	Matrix Ed();
	Matrix operator = (const Matrix matrix);
	Matrix operator + ( Matrix& matr);
	Matrix operator - ( Matrix& matrix);
	Matrix operator * ( Matrix& matrix);
	Matrix operator * (const double number);
    bool operator == (const Matrix& matrix);
	Matrix transposed_matrix();
	void swap_str(int k, int r);
	int return_n();
	double return_value_of_matrix(int str, int col);
	void recording_value(int str, int col, double value);
	void get_matrix(std::istream& sth);
	void print_matrix();
	void create_memory(int n, int m);
	void delete_memory();
	int return_m();
};




class Vector {

private:
	int size;
	double* vec;
public:
	Vector(int size = -1, std::istream& sth = std::cin);
	double Norma_1();
	double Norma_2();
	double Norma_8();
	void print();
	Vector(const Vector& vector);
	~Vector();
	Vector operator * (Matrix& matrix);
	bool operator == ( Vector& v);
	Vector operator + ( Vector& v);
	Vector operator - ( Vector& v);
	double operator * ( Vector& v);
	Vector operator = (const Vector vector);
	Vector operator * ( double number);
	void get_size(std::istream& sth);
	void get_vector(std::istream& sth);
	void swap_vector(int k, int r);
	void change_value_vector(int index_current, int index_dif, double mat_dif, double mat_current);
	void value_for_solve(int index_current, double mat_current);
	void solve_for_vector(int  index_current, double mat_current);
	int return_size();
	double return_values();
	double return_one_value(int i);
	void recording_value(int i, double value);
	void create_memory(int size);
	void delete_memory();
};