#include "VectorsOperation.h"
#include <algorithm>

std::vector<double> SubtractionVector(std::vector<double> a, std::vector<double> b)
{
	std::vector<double> res;
	for (int i = 0; i < a.size(); ++i)
		res.push_back(a[i] - b[i]);

	return res;
}

double FirstVectorNorm(std::vector<double> a)
{
	return abs(*std::max_element(a.begin(), a.end(), [](double x, double y) {return abs(x) < abs(y); }));
}

double SecondVectorNorm(std::vector<double> a)
{
	double sum = 0;

	for (int i = 0; i < a.size(); ++i)
		sum += std::abs(a[i]);

	return sum;
}

double ScalarMult(std::vector<double> a, std::vector<double> b)
{
	double sum = 0;

	for (int i = 0; i < a.size(); ++i)
		sum += a[i] * b[i];

	return sum;
}

double ThirdVectorNorm(std::vector<double> a)
{
	return sqrt(ScalarMult(a, a));
}

std::vector<double> AdditionVector(std::vector<double> A, std::vector<double> B) // ÿ íå ïîìíþ åñòü ëè ýòî â òâîåé áèáå èëè íåò
{
	std::vector<double> res;
	for (int i = 0; i < A.size(); ++i)
		res.push_back(A[i] + B[i]);
	return res;
}

std::vector<double> MultVectorNum(std::vector<double> A, double B) // ÿ íå ïîìíþ åñòü ëè ýòî â òâîåé áèáå èëè íåò
{
	std::vector<double> res;
	for (int i = 0; i < A.size(); ++i)
		res.push_back(A[i] * B);
	return res;
}