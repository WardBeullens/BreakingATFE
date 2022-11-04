#ifndef FORM_H
#define FORM_H

#include "matrix.h"
#include "minrank.h"

typedef std::array<long long, FORM_LEN> Form;

void randomize_form(Form &f);

Form act_on_form(const Form &f, const Matrix &S);

Matrix form_to_mat(const Form f, const Vector &v);
void brute_force_rank(const Form f, const unsigned int r, Vector &v);

std::vector<Vector> enumerate_low_rank(const Form &f);
std::vector<Vector> enumerate_low_rank_neighbours(const Form &f, const Vector &v);
Vector find_rank_4(const Form &f);
Vector F(const Form &f, const Vector v);
std::vector<Vector> compute_F_chain(const Form &f, const Vector &start, size_t Len);
Vector extract_signature(const Form &f, const std::vector<Vector> &chain, size_t pos);
std::pair<Vector,Vector> find_colliding_pair(const Form &f1, const Form &f2);
Matrix find_isomorphsism(const Form &f1, const Vector &v1, const Form &f2, const Vector &v2);
Vector random_neighbour(const Form &f, const Vector &v);

#endif 