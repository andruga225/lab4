// lab4.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include "VectorsOperation.h"
#include "MatrixOperation.h"

using namespace std;

int variant=4;
ofstream fout("output.txt");

vector<vector<double>> Jacobian(double x, double y)
{
	
	vector<vector<double>> Jac(2, (vector<double>(2)));
	if (variant == 4) {
		Jac[0][0] = sin(x);
		Jac[0][1] = 0;
		Jac[1][0] = 0;
		Jac[1][1] = 0.5 * cos(y - 0.5);
	}else
	{
		Jac[0][0] = 0;
		Jac[0][1] = -cos(y - 1);
		Jac[1][0] = cos(x + 1);
		Jac[1][1] = 0;
	}

	return Jac;
}

double F1(double x, double y)
{
	if (variant == 4)
		return cos(x) + y - 1.5;
	else
		return x + sin(y - 1) - 1.3;
}

double F2(double x, double y)
{
	if (variant == 4)
		return 2 * x - sin(y - 0.5) - 1;
	else
		return y - sin(x + 1) - 0.8;
}

double FF(double x, double y)
{
	return pow(F1(x,y),2) + pow(F2(x,y),2);
}

double DFFX(double x, double y)
{
	if (variant == 4)
		return 2 * (cos(x) + y - 1.5) * (-sin(x)) + 2 * (2 * x - sin(y - 0.5) - 1) * 2;
	else
		return 2 * sin(y - 1) + 2 * x - 2.6 - 2 * (y - 0.8) * cos(x + 1) + sin(2 * x + 2);
}

double DFFY(double x, double y)
{
	if (variant == 4)
		return 2 * (cos(x) + y - 1.5) + 2 * (2 * x - sin(y - 0.5) - 1) * (-cos(y - 0.5));
	else
		return sin(2 * y - 2) + 2 * cos(y - 1) * (x - 1.3) + 2 * y - 1.6 - 2 * sin(x + 1);
}

void SimpleIteration(double x, double y)
{
	vector<vector<double>> Jac = Jacobian(x, y);
	int itr = 0;
	double error=1, nevNorm = 1;

	do
	{
		itr++;

		double xNew, yNew;
		if (variant == 4) {
			xNew = 0.5+0.5*sin(y-0.5);
			yNew = 1.5-cos(x);
		}else
		{
			xNew = 1.3 - sin(y - 1);
			yNew = 0.8 + sin(x + 1);

		}

		Jac = Jacobian(xNew, yNew);
		double q = cubic_norm(Jac);
		vector<double> nev;

		nev.push_back(abs(F1(xNew,yNew)));
		nev.push_back(abs(F2(xNew,yNew)));
		nevNorm = ThirdVectorNorm(nev);

		error = q * ThirdVectorNorm({ xNew - x,yNew - y })/(1-q);
		
		fout <<fixed<<setw(2)<< itr <<" | "<<fixed<<right<<setw(10)<<setprecision(7)<<xNew<<" | "<<yNew<< " || " << fixed<<setw(20)<<setprecision(15)<<nevNorm << " | " << F1(xNew,yNew)<<" | "<<F2(xNew,yNew)<<" | " << fixed << setw(20) << setprecision(15) << q << endl;

		x = xNew;
		y = yNew;

	}while (nevNorm > 0.0001 && error > 0.0001);

}

vector<vector<double>> dF(double x, double y)
{
	vector<vector<double>> f(2, (vector<double>(2)));

	if (variant == 4) {
		f[0][0] = -sin(x);
		f[0][1] = 1;
		f[1][0] = 2;
		f[1][1] = -cos(y-0.5);
	}else
	{
		f[0][0] = 1;
		f[0][1] = cos(y - 1);
		f[1][0] = -cos(x + 1);
		f[1][1] = 1;
	}

	return f;
}

void NewtonMethod(double x, double y)
{
	vector<vector<double>> Jac = Jacobian(x, y);
	int itr = 0;
	double error = 1, nevNorm = 1;
	vector<vector<double>> f = dF(x, y);

	do
	{
		itr++;

		f = dF(x, y);
		f = ReverseMatrix(f);
		double f1 = F1(x, y);
		double f2 = F2(x, y);

		double xNew = x - (f[0][0] * f1 + f[0][1] * f2);
		double yNew = y - (f[1][0] * f1 + f[1][1] * f2);

		Jac = Jacobian(xNew, yNew);
		double q = cubic_norm(Jac);

		vector<double> nev;

		nev.push_back(abs(F1(xNew, yNew)));
		nev.push_back(abs(F2(xNew, yNew)));
		nevNorm = ThirdVectorNorm(nev);

		error = q * ThirdVectorNorm({ xNew - x,yNew - y }) / (1 - q);

		fout << fixed << setw(2) << itr << " | " << fixed << right << setw(10) << setprecision(7) << xNew << " | " << yNew << " || " << fixed << setw(20) << setprecision(15) << nevNorm << " | " << F1(xNew, yNew) << " | " << F2(xNew, yNew) << " | " << fixed << setw(20) << setprecision(15) << q << endl;

		x = xNew;
		y = yNew;

	} while (nevNorm > 0.0001 && error > 0.0001||itr<3);
}


void Gradient(double x, double y)
{

	vector<vector<double>> Jac = Jacobian(x, y);
	int itr = 1;
	double error = 1, nevNorm = 1;
	double xNew, yNew;

	do
	{
		double alpha = 1;
		double lambda = 0.5;

		double dx = DFFX(x, y);
		double dy = DFFY(x, y);

		while (FF(x - alpha * dx, y - alpha * dy) >= FF(x, y))
			alpha *= lambda;

		xNew = x - alpha * dx;
		yNew = y - alpha * dy;

		Jac = Jacobian(xNew, yNew);
		double q = cubic_norm(Jac);

		vector<double> nev;

		nev.push_back(abs(F1(xNew, yNew)));
		nev.push_back(abs(F2(xNew, yNew)));

		nevNorm = ThirdVectorNorm(nev);
		error = ThirdVectorNorm({ xNew - x,yNew - y }) * q / (1 - q);

		fout <<fixed<< setw(2)<< itr << " | " <<fixed<<setw(10)<<setprecision(7) <<xNew << " | " << yNew<<" || "<<alpha<<" | "<<fixed<<setw(20)<<setprecision(15)<<nevNorm<<" | "<< F1(xNew, yNew) << " | " << F2(xNew, yNew) <<" | "<<nevNorm<< endl;
		
		itr++;
		x = xNew;
		y = yNew;

	} while (nevNorm > 0.0001 && error > 0.0001);
}

int main()
{
	if (variant == 4)
		fout << "x0 = " << 0.5 << " y0 = " << 0.5 << endl;
	else
		fout << "x0 = " << 0.6 << " y0 = " << 1.8 << endl;

	fout << "Метод простой итерации" << endl;
	if (variant == 4) {
		fout << "Fi1(x,y)=1.5-cos(x)" << endl;
		fout << "Fi2(x,y)=0.5+0.5sin(y-0.5)" << endl;
	}
	else
	{
		fout << "Fi1(x,y)=1.3 - sin(y-1)" << endl;
		fout << "Fi2(x,y)=0.8 + sin(x+1)" << endl;
	}

	fout << "Якобиан" << endl;

	if(variant==4)
	{
		fout << "sin(x) 0.0" << endl;
		fout << "0.0 -cos(y-0.5)" << endl;
	}
	else
	{
		fout << "0.0 -cos(y-1)" << endl;
		fout << "cos(x+1) -0.0" << endl;
	}

	fout << "Значение" << endl;

	vector<vector<double>> Jac;

	if (variant == 4)
		Jac = Jacobian(0.5, 0.5);
	else
		Jac = Jacobian(0.6, 1.8);

	for(int i=0;i<Jac.size();++i)
	{
		for (int j = 0; j < Jac[i].size(); ++j)
			fout << fixed << setprecision(4) << Jac[i][j] << " ";
		fout << endl;
	}

	fout << "Норма" << endl;

	fout << fixed << setprecision << cubic_norm(Jac) << endl;
	fout << "Itr x y Норма невязки F1 F2 Норма якобиана" << endl;
	if (variant == 4)
		SimpleIteration(0.5, 0.5);
	else
		SimpleIteration(0.6, 1.8);

	fout << "Метод Ньютона" << endl;
	fout << "Матрица производных" << endl;

	if(variant==4)
	{
		fout << "-sin(x) 1" << endl;
		fout << "2 -cos(y-0.5)" << endl;
	}
	else
	{
		fout << "1 cos(y-1)" << endl;
		fout << "-cos(x+1) 1" << endl;
	}

	fout << "Itr x y Норма невязки F1 F2" << endl;
	if (variant == 4)
		NewtonMethod(0.5, 0.5);
	else
		NewtonMethod(0.6, 1.8);


	fout << "Метод градиентного спуска" << endl;
	fout << "Itr x y Alfa Норма невязки F1 F2 FF k" << endl;
	if (variant == 4)
		Gradient(0.5, 0.5);
	else
		Gradient(0.6, 1.8);
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
