#pragma once
#pragma once
#include <vector>

const double eps = 0.00000000001;
std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>>);
std::vector<std::vector<double>> MultMatrix(std::vector<std::vector<double>>, std::vector<std::vector<double>>);
std::vector<double> MultMatrixVector(std::vector<std::vector<double>>, std::vector<double>);
double EuclideanNorm(std::vector<std::vector<double>>);
std::vector<std::vector<double>> SubtractionMatrix(std::vector<std::vector<double>>, std::vector<std::vector<double>>);
double cubic_norm(std::vector<std::vector<double>>);
double octahedral_norm(std::vector<std::vector<double>>);
std::vector<std::vector<double>> ReverseMatrix(std::vector<std::vector<double>>);