#ifndef MINRANK_H
#define MINRANK_H

#include "matrix.h"

std::vector<std::vector<long long>> enumerate_low_rank(const std::vector<Matrix> input_matrices);
int salvage_attempt(Matrix &CM, const std::vector<Matrix> Matrices, std::vector<Vector> &sols);

#endif