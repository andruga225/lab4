// lab4.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <iomanip>
#include "VectorsOperation.h"
#include "MatrixOperation.h"

using namespace std;

vector<vector<double>> Jacobian(double x, double y)
{
	//Тут для его варика пока что
	vector<vector<double>> Jac(2, (vector<double>(2)));
	Jac[0][0] = 0;
	Jac[0][1] = cos(y+0.5);
	Jac[1][0] = sin(x-2);
	Jac[1][1] = 0;

	return Jac;
}

double F1(double x, double y)
{
	return x - sin(0.5 + y) + 1;
}

double F2(double x, double y)
{
	return y + cos(x-2);
}

double phix(double y)
{
	return sin(0.5 + y) - 1;
}
double phiy(double x)
{
	return -cos(x - 2);
}
void SimpleIteration(double x, double y)
{
	vector<vector<double>> Jac = Jacobian(x, y);
	int itr = 0;
	double error=1, nevNorm = 1;

	do
	{
		itr++;

		double xNew = -1 + sin(y+0.5);
		double yNew = -cos(x-2);

		Jac = Jacobian(xNew, yNew);
		double q = cubic_norm(Jac);
		vector<double> nev;

		nev.push_back(abs(F1(xNew,yNew)));
		nev.push_back(abs(F2(xNew,yNew)));
		nevNorm = ThirdVectorNorm(nev);

		error = q * ThirdVectorNorm({ xNew - x,yNew - y })/(1-q);
		
		cout << itr <<" " <<nevNorm<<" "<<fixed<<setprecision(15)<<q<< endl;

		x = xNew;
		y = yNew;

	}while (nevNorm > 0.0001 && error > 0.0001);

}
int main()
{
    std::cout << "Hello World!\n";
	SimpleIteration(-0.1, 0.5);
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
