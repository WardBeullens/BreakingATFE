
#ifndef MATRIX_H
#define MATRIX_H

#include <cassert>
#include <iostream>
#include "stdint.h"
#include "cstdlib"
#include <time.h>  
#include <array>
#include <vector>
#include <unordered_set>
#include <stack>

#ifdef _MSC_VER
# include <intrin.h>
#else
# include <x86intrin.h>
#endif

#define TIC uint64_t TIC_TOC_start = __rdtsc();
#define TOC(A) std::cout << #A << ": " << __rdtsc()-TIC_TOC_start << std::endl; TIC_TOC_start = __rdtsc();

#define N 9
#define FORM_LEN (N*(N-1)*(N-2)/6)
#define Q 524287 //31 // 101 // 257 // 997 // 2003 // 10007 // 50021 // 524287
#define RANK (N-5)
#define N_CHOOSE_R (N*(N-1)*(N-2)*(N-3)*(N-4)/120)
#define N_CHOOSE_R_PLUS_ONE (N*(N-1)*(N-2)*(N-3)/24)

#define Data(i,j) data[((i)*col_length) + (j)]
#define MData(mat,i,j) mat.data[(((i)*(mat.col_length))) + (j)]

using Vector = std::vector<long long>;

long long inverse_mod_Q(long long a);
std::ostream& operator<<(std::ostream& os, const Vector &vec);

template<size_t LEN> void print_array(const std::array<long long, LEN> array);
template<size_t LEN> int normalize_array(std::array<long long, LEN> &array);

void randomize_vector(Vector &vec);
int normalize_vector(Vector &v);


long long operator*(const Vector &v1, const Vector &v2);
Vector operator*(long long c, const Vector &v);

Vector& operator+=(Vector &v_acc, const Vector &v);
Vector operator+(const Vector &v1, const Vector &v2);


class Matrix
{
private:
public:
    size_t rows;
    size_t cols;
    size_t col_length;
    Vector data;

    Matrix(size_t rows = 0, size_t cols = 0);

    void reduce();
    void randomize();
    void zero();
    void identity();
    void EF();
    void RREF();
    size_t rank() const;
    size_t rank_EF() const;

    const int pivot_row(const size_t col);
    const std::vector<Vector> kernel();
    const Vector random_kernel_vector_EF();
    const Vector kernel_vector_EF(size_t i = 0);
    const std::vector<Vector> kernel_EF();
    Vector row(size_t r) const;
    Vector col(size_t c) const;
    void set_row(size_t i, const Vector &v);
    void set_col(size_t i, const Vector &v);

    Vector operator*(const Vector&) const;

    Matrix& operator+=(const Matrix& rhs);
    Matrix& operator-=(const Matrix& rhs);
    friend std::ostream& operator<<(std::ostream& os, const Matrix &Mat);
    friend void add_mul(long long scalar, const Matrix& summand, Matrix& accumulator);
};

Matrix random_invertible_matrix(size_t size);
int multiply_by_inverse(const Matrix &mat, const Matrix &in, Matrix &out);

#endif