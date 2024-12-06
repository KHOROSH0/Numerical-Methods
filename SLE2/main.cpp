#include "Matrix_Vector.h"
#include "MPI_ZEIDEL_LU_QR.h"
#include "Manip.h"
#include "istream"
#include "ostream"
#include "fstream"
#include <Vector>
#include <string>
using namespace std;
int main() {
	std::vector <Matrix> A;
	std::vector <Vector> b;
	std::string name_of_file = "Matrix_Vector.txt";
	std::ifstream fin;
	fin.open(name_of_file);
	int n = 0, m = 0, size = -1;
	for (int i = 0; i < 5; i++) {
		Matrix A_{ 0,0, fin };
		Vector B_{ -1, fin };
		A.push_back(A_);
		b.push_back(B_);
	}
	fin.close();
	std::vector <Matrix> A_result;
	std::vector <Vector> b_result;
	std::vector <Vector> bb_result;
	std::string name_of_file2 = "Matrix_Vector_results.txt";
	fin.open(name_of_file2);
	int i = 0;
	while (fin.peek() != EOF) {
		if (i < 5) {
			Vector B_{ -1, fin };
			b_result.push_back(B_);
		}
		else {
			Vector B_{ -1, fin };
			bb_result.push_back(B_);
		}
		i++;
	}
	int r = 0;
	std::vector <Matrix> AA;
	std::vector <Vector> bb;
	std::vector <double> eps1{ 1e-3, 1e-6 };
	for (int i = 0; i < 10; i++) {
		for (int l = 0; l < 2; l++) {
			Matrix A_{ i + 4,i + 4 };
			Vector B_{ i + 4 };
			for (int j = 0; j < A_.return_m(); j++) {
				if (j != A_.return_m() - 1) {
					B_.recording_value(j, -1);
				}
				else {
					B_.recording_value(j, 1);
				}
				for (int k = 0; k < A_.return_n(); k++) {
					if (j == k) {
						A_.recording_value(j, k, 1 + eps1[l] * 13);
					}
					else if (k < j) {
						A_.recording_value(j, k, 0 + eps1[l] * 13);
					}
					else {
						A_.recording_value(j, k, -1 - eps1[l] * 13);
					}
				}
			}
			AA.push_back(A_);
			bb.push_back(B_);
			r = i + l;
		}
	}
	vector <double> eps{ 1e-3, 1e-4, 1e-5, 1e-6 };
	vector <string> head1 = { " number of test", " x_from_lib", "eps", "MPI", "x", "delta                      ", "k", "Zidel","x", "delta", "k", "LU","x", "delta", "QR", "x", "delta" };
	Writer writer(head1);
	for (int i(0); i < 4; i++) {
		int test = 1 + i;
		for (int j(0); j < 4; j++) {
			vector <string> value = {
			 to_string(test), toString(b_result[i]), to_string(eps[j]),
			 " ", toString(MPI(A[i], b[i], eps[j])), to_string((MPI(A[i], b[i], eps[j]) - b_result[i]).Norma_2()), to_string(iter_of_MPI(A[i], b[i], eps[j])),
			 " ", toString(Zeidel(A[i], b[i], eps[j])), to_string((Zeidel(A[i], b[i],eps[j]) - b_result[i]).Norma_2()), to_string(iter_of_Zeidel(A[i], b[i], eps[j])),
			 " ", toString(LU(A[i], b[i])), to_string((LU(A[i], b[i]) - b_result[i]).Norma_2()),
			 " ", toString(QR(A[i], b[i])), to_string((QR(A[i], b[i]) - b_result[i]).Norma_2())
			};
			writer.add(value);
		}
	}
	writer.print();
	vector <string> head2 = { " number of test"," dim","eps2", " x_from_lib", "eps", "MPI", "x", "delta", "k", "Zidel","x", "delta", "k", "LU","x", "delta", "QR", "x", "delta" };
	Writer writer1(head2);
	int dim = 4;
	for (int i = 0; i < 18; i++) {
		int test = 5;
		for (int k = 0; k < 2; k++) {
			for (int j = 0; j < 4; j++) {
				vector <string> value = {
					to_string(test),to_string(dim),to_string(eps1[k]), toString(bb_result[i]), to_string(eps[j]),
					" ", toString(MPI(AA[i], bb[i], eps[j])), to_string((MPI(AA[i], bb[i], eps[j]) - bb_result[i]).Norma_2()), to_string(iter_of_MPI(AA[i], bb[i], eps[j])),
					" ", toString(Zeidel(AA[i], bb[i], eps[j])), to_string((Zeidel(AA[i], bb[i],eps[j]) - bb_result[i]).Norma_2()), to_string(iter_of_Zeidel(AA[i], bb[i], eps[j])),
					" ", toString(LU(AA[i], bb[i])), to_string((LU(AA[i], bb[i]) - bb_result[i]).Norma_2()),
					" ", toString(QR(AA[i], bb[i])), to_string((QR(AA[i], bb[i]) - bb_result[i]).Norma_2())
				};
				writer1.add(value);
			}
		}
		dim++;
	}
	writer1.print();
}