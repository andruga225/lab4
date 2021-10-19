#include "MatrixOperation.h"

std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>> a)
{
    std::vector<std::vector<double>> res(a.size(), (std::vector<double>(a[0].size())));

    for (int i = 0; i < a.size(); ++i)
        for (int j = 0; j < a[i].size(); ++j)
            res[j][i] = a[i][j];

    return res;
}

std::vector<std::vector<double>> MultMatrix(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b)
{
    std::vector<std::vector<double>> res(a.size(), (std::vector<double>(a[0].size())));

    for (int i = 0; i < a.size(); ++i)
        for (int j = 0; j < a.size(); ++j)
        {
            double s = 0;

            for (int k = 0; k < a.size(); ++k)
                s += a[i][k] * b[k][j];

            res[i][j] = s;
        }

    return res;
}

std::vector<double> MultMatrixVector(std::vector<std::vector<double>> a, std::vector<double> b)
{
    std::vector<double> result(b.size());

    for (int i = 0; i < a.size(); ++i)
    {
        double s = 0;
        for (int j = 0; j < b.size(); ++j)
            s += a[i][j] * b[j];

        result[i] = s;
    }

    return result;
}

void search_max(int& i, int& j, std::vector<std::vector<double>> A)
{
    double max = LONG_MIN;

    for (int x = 0; x < A.size(); x++)
    {
        for (int y = 0; y < A[x].size(); y++)
        {
            if (abs(A[x][y]) > max && x != y)
            {
                max = abs(A[x][y]);
                i = x;
                j = y;
            }
        }
    }
}

double summ(std::vector<std::vector<double>> A)
{
    double sum = 0;
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            if (i != j)
            {
                sum += pow(A[i][j], 2.0);
            }
        }
    }
    return sum;
}

double EuclideanNorm(std::vector<std::vector<double>> A)
{
    int x, y;

    std::vector<std::vector<double>> T(A.size(), std::vector<double>(A.size()));
    std::vector<std::vector<double>> TT;

    std::vector<std::vector<double>> AA;


    for (int i = 0; i < T.size(); i++)
    {
        for (int j = 0; j < T.size(); j++)
        {
            T[i][j] = 0;
            if (i == j) T[i][j] = 1;
        }
    }

    double alpha;

    while (abs(summ(A)) > eps) {

        for (int i = 0; i < T.size(); i++)
        {
            for (int j = 0; j < T.size(); j++)
            {
                T[i][j] = 0;
                if (i == j) T[i][j] = 1;
            }
        }

        search_max(x, y, A);
        double d = A[x][y];

        if (abs(A[x][x] - A[y][y]) < eps) {
            alpha = atan(1);
        }
        else
        {
            alpha = abs(atan(2 * A[x][y] / (A[x][x] - A[y][y]))) / 2;
        }

        T[x][x] = cos(alpha);
        T[x][y] = -sin(alpha);
        T[y][x] = sin(alpha);
        T[y][y] = cos(alpha);

        TT = Transpose(T);

        AA = MultMatrix(TT, A);

        A = MultMatrix(AA, T);
    }

    double min = LONG_MAX;
    double max = LONG_MIN;

    for (int i = 0; i < A.size(); i++)
    {
        if (A[i][i] > max)
            max = A[i][i];

        if (A[i][i] < min)
            min = A[i][i];
    }

    return sqrt(max);
}

std::vector<std::vector<double>> SubtractionMatrix(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b)
{
    std::vector<std::vector<double>> res(a.size(), (std::vector<double>(a.size())));

    for (int i = 0; i < res.size(); ++i)
        for (int j = 0; j < res.size(); ++j)
            res[i][j] = a[i][j] - b[i][j];

    return res;
}

double cubic_norm(std::vector<std::vector<double>> mass)
{
    double max = LONG_MIN;
    double summ = 0;

    for (int i = 0; i < mass.size(); i++)
    {
        for (int j = 0; j < mass[i].size(); j++)
        {
            summ += abs(mass[i][j]);
        }

        if (summ - max > eps)
            max = summ;

        summ = 0;
    }

    return max;
}

double octahedral_norm(std::vector<std::vector<double>> mass)
{
    double max = LONG_MIN;
    double summ = 0;

    for (int i = 0; i < mass[0].size(); i++)
    {
        for (int j = 0; j < mass.size(); j++)
        {
            summ += abs(mass[j][i]);
        }
        if (summ - max > eps) max = summ;
        summ = 0;
    }
    return max;
}

std::vector<std::vector<double>> ReverseMatrix(std::vector<std::vector<double>> a)
{
    std::vector<std::vector<double>> e(a.size(), (std::vector<double>(a.size())));
    std::vector<std::vector<double>> RevA = a;

    for (int i = 0; i < e.size(); ++i)
        e[i][i] = 1;

    for(int k=0;k<a.size();++k)
    {
        double temp = RevA[k][k];

        for(int j=0;j<a.size();++j)
        {
            RevA[k][j] /= temp;
            e[k][j] /= temp;
        }

        for(int i=k+1;i<a.size();++i)
        {
            temp = RevA[i][k];

            for(int j=0;j<a.size();++j)
            {
                RevA[i][j] -= RevA[k][j] * temp;
                e[i][j] -= e[k][j] * temp;
            }
        }
    }

    for(int k=a.size()-1;k>0;--k)
    {
	    for(int i=k-1;i>-1;--i)
	    {
            double temp = RevA[i][k];
            for(int j=0;j<a.size();++j)
            {
                RevA[i][j] -= RevA[k][j] * temp;
                e[i][j] -= e[k][j] * temp;
            }
	    }
    }

    return e;
}

