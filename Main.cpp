#include <stdio.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;

const double ZERO = 1E-20; // ?????, ?? ?????? ?????? ZERO, ??????? ?????
const int maxiter = 100; // ???????? ???????? ?????? ????????/????????? ???????
int funcComp = 0; // ???-?? ?????????? ???????

double r = 10; // ??????????? ?????? ? ??? ????????? ????????
vector<double> xk = { -10, -10 }; // ????????? ???????????
double epsF = 1e-7; // ???????? ?????? ????????/????????? ???????
const double M = 1E+7; // ??????? ??????????????? ? epsF
int flag = 1; // 1 ??? ?????????, 0 ??? ????????

//-----------------------------------------------------------------------------------------------
// ???????, ??????????????? ???????, ????????? ?????????
//-----------------------------------------------------------------------------------------------

double change_r(double r)
{
	return r / 1000;
}

// ???????? ???????
double func(double x, double y)
{
	funcComp++;
	return 4 * (y - x) * (y - x) + 3 * (x - 1) * (x - 1);
}

// ?????? ? (?????????)
double func_h(double x, double y)
{
	return x + 1 - y;
}

// ?????? ? (???????????)
double func_g(double x, double y)
{
	return x + y + 1;
}

// ???????? ??????? ??? ?????? ?
double func_H(double x, double y)
{
	//return abs(func_h(x, y));
	return func_h(x, y) * func_h(x, y);
	//return pow(abs(func_h(x, y)), 5);
}

// ???????? ??????? ??? ?????? ?
double func_G(double x, double y)
{
	return 0.5 * (func_g(x, y) + abs(func_g(x, y)));
	//return 0.5 * (func_g(x, y) + abs(func_g(x, y))) * 0.5 * (func_g(x, y) + abs(func_g(x, y)));
	//return pow(0.5 * (func_g(x, y) + abs(func_g(x, y))), 5);
}

// ????????? ???????
double func_Gb(double x, double y)
{
	//return -1. / func_g(x, y);
	return -log(-func_g(x, y));
}

// ??????????????? ???????
double GetQ(vector<double>& xk)
{
	// ?????? ?
	//return func(xk[0], xk[1]) + r * (func_G(xk[0], xk[1]));
	// ?????? ?
	return func(xk[0], xk[1]) + r * (func_H(xk[0], xk[1]));
}

double GetFine(vector<double>& xk)
{
	// ?????? ?
	//return 1 * (func_G(xk[0], xk[1]));
	// ?????? ?
	return 1 * (func_H(xk[0], xk[1]) );
}

double GetQb(vector<double>& xk)
{
	return func(xk[0], xk[1]) + r * (func_Gb(xk[0], xk[1]));
}

double GetBarrier(vector<double>& xk)
{
	return 1 * (func_Gb(xk[0], xk[1]));
}

//-----------------------------------------------------------------------------------------------
// ????? ???????? ??????? ?????????????? ????????????? ?? ????????
//-----------------------------------------------------------------------------------------------
const int N = 2,
MAXITER = 100;
const double EPS = 1E-10;
vector<vector<double>> D; // ?????????? ?????? ?????????
vector<double> xMiddle, xNew;
double fMin, fMax;
int xMin, xMax;
double alpha = 1,
gamma = 2.5,
beta = 0.5;

int Prepare(vector<double>& x0)
{
	D.resize(N + 1);
	xMiddle.resize(N);
	xMiddle = { 0, 0 };
	xNew.resize(N);
	xNew = { 0, 0 };

	double t = 5; // ?????????? ????? ?????????
	double d1 = (sqrt(N + 1) + N - 1) * t / (N * sqrt(2)),
		d2 = t / ((sqrt(N + 1) - 1) * N * sqrt(2));
	D[0].resize(N);
	D[0] = x0;
	for (int i = 1; i < N + 1; i++)
	{
		D[i].resize(N);
		for (int j = 0; j < N; j++)
		{
			D[i][j] = D[0][j] + d2;
		}
		D[i][i - 1] = D[0][i - 1] + d1;
	}

	double koef = 1;
	while (flag && (func_g(D[1][0], D[1][1]) > 0 || func_g(D[2][0], D[2][1]) > 0))
	{
		koef /= 2;
		for (int i = 0; i < N; i++)
		{
			D[1][i] = D[0][i] + koef * d2;
			D[2][i] = D[0][i] + koef * d2;
		}
	}

	return 0;

}

int FindMiddle()
{
	for (int j = 0; j < N; j++)
	{
		double sum = 0;
		for (int i = 0; i < N + 1; i++)
		{
			sum += D[i][j];
		}
		sum -= D[xMax][j];
		xMiddle[j] = sum / N;
	}
	return 0;
}

double FindMinMax()
{
	double ftemp = flag ? GetQb(D[0]) : GetQ(D[0]),
		fmiddle = flag ? GetQb(xMiddle) : GetQ(xMiddle);
	fMax = ftemp;
	fMin = ftemp;
	xMax = 0;
	xMin = 0;
	double sum = 0;
	for (int i = 1; i < N + 1; i++)
	{
		ftemp = flag ? GetQb(D[i]) : GetQ(D[i]);
		sum += (ftemp - fmiddle) * (ftemp - fmiddle);
		if (ftemp > fMax)
		{
			fMax = ftemp;
			xMax = i;
		}
		else if (ftemp < fMin)
		{
			fMin = ftemp;
			xMin = i;
		}
	}
	return sum;
}


int Reflection() // ?????????
{
	for (int i = 0; i < N; i++)
	{
		xNew[i] = xMiddle[i] + alpha * (xMiddle[i] - D[xMax][i]);
	}

	double koef = alpha;
	while (flag && func_g(xNew[0], xNew[1]) > 0)
	{
		koef /= 2;
		for (int i = 0; i < N; i++)
		{
			xNew[i] = xMiddle[i] + koef * (xMiddle[i] - D[xMax][i]);
		}
	} 
	return 0;
}

int Expansion(double fTemp) // ??????????
{
	vector<double> xTemp(N);
	for (int i = 0; i < N; i++)
	{
		xTemp[i] = xMiddle[i] + gamma * (xNew[i] - xMiddle[i]);
	}

	if (flag && func_g(xTemp[0], xTemp[1]) > 0)
	{
		D[xMax] = xNew;
	}
	else
	{
		double Q_temp = flag ? GetQb(xTemp) : GetQ(xTemp);
		if (Q_temp < fTemp)
		{
			D[xMax] = xTemp;
		}
		else
		{
			D[xMax] = xNew;
		}
	}
	return 0;
}

int Contraction(double fTemp) // ??????
{
	vector<double> xTemp(N);
	for (int i = 0; i < N; i++)
	{
		xTemp[i] = xMiddle[i] + beta * (D[xMax][i] - xMiddle[i]);
	}

	double koef = beta;
	while (flag && func_g(xTemp[0], xTemp[1]) > 0)
	{
		koef /= 2;
		for (int i = 0; i < N; i++)
		{
			xTemp[i] = xMiddle[i] + beta * (D[xMax][i] - xMiddle[i]);
		}
	}

	double fTemp2 = flag ? GetQb(xTemp) : GetQ(xTemp);
	if (fTemp2 < fTemp)
		D[xMax] = xTemp;
	else
		D[xMax] = xNew;
	return 0;

}

int Shrink() // ????????
{
	for (int i = 0; i < N + 1; i++)
	{
		if (i != xMin)
		{
			for (int j = 0; j < N; j++)
			{
				D[i][j] = D[xMin][j] + 0.5 * (D[i][j] - D[xMin][j]);
			}
		}
	}
	return 0;

}

int Algorithm_Polytype(vector<double>& x0)
{
	Prepare(x0);
	int iter;
	double stop = 10;
	FindMinMax();

	for (iter = 1; iter < MAXITER && stop > EPS; iter++)
	{
		FindMiddle();

		Reflection();
		double fTemp = flag ? GetQb(xNew) : GetQ(xNew);
		if (fTemp > fMax)
		{
			Shrink();
		}
		else if (fTemp < fMax)
		{
			if (fTemp < fMin)
			{
				Expansion(fTemp);
			}
			else
			{
				Contraction(fTemp);
			}
		}

		stop = FindMinMax();
		stop = sqrt(stop);
		stop /= (N + 1);

	}

	x0 = D[xMin];
	return iter;
}

//-----------------------------------------------------------------------------------------------
// ????? ???????? ???????
//-----------------------------------------------------------------------------------------------
void Fine()
{
	double fine = GetFine(xk);
	ofstream fout;
	fout.open("out.txt");
	fout.imbue(locale("Russian"));
	fout << setprecision(15);
	int ii;
	cout << "fine=" << fine << " xk: " << xk[0] << " " << xk[1] << endl;
	fout << "fine=" << fine << " xk: " << xk[0] << " " << xk[1] << endl;
	for (ii = 0; ii < maxiter && abs(fine) > epsF; ii++)
	{
		r = change_r(r);
		Algorithm_Polytype(xk);
		fine = r * GetFine(xk);
		cout << "fine=" << fine << " xk: " << xk[0] << "\t" << xk[1] << "\t" << r << "\t" << xk[0] + xk[1] << "\t" << func(xk[0], xk[1]) << endl;

		fout << "fine=" << fine << " xk: " << xk[0] << " " << xk[1] << endl;
	}

	fout << "?????: \n" << ii << '\t' << funcComp << '\t' << func(xk[0], xk[1]) << '\t' << xk[0] << '\t' << xk[1];
	cout << "f = " << func(xk[0], xk[1]) << endl;
	cout << "Fine Iter: " << ii << endl;
	// ???-?? ?????????? ???????
	cout << "func compute: " << funcComp << endl;
}

//-----------------------------------------------------------------------------------------------
// ????? ????????? ???????
//-----------------------------------------------------------------------------------------------
void Barrier()
{
	ofstream fout;
	fout.open("out.txt");
	fout.imbue(locale("Russian"));
	fout << setprecision(15);
	double barrier = 1;

	int ii;
	for (ii = 0; ii < maxiter && abs(barrier) > epsF; ii++)
	{
		r = change_r(r);		
		Algorithm_Polytype(xk); //????? ????? ???????? xk

		if (abs(func_g(xk[0], xk[1])) < ZERO)
		{
			barrier = r / M;
		}
		else
		{
			barrier = r * GetBarrier(xk);
		}
		cout << "bar=" << barrier << " xk: " << xk[0] << "\t" << xk[1] << "\t" << r << "\t" << xk[0] + xk[1] << "\t" << func(xk[0], xk[1]) << endl;
		fout << "bar=" << barrier << " xk: " << xk[0] << "\t" << xk[1] << endl;
	}
	fout << "?????: \n" << ii << '\t' << funcComp << '\t' << func(xk[0], xk[1]) << '\t' << xk[0] << '\t' << xk[1];
	cout << "f = " << func(xk[0], xk[1]) << endl;
	cout << "Barrier Iter: " << ii << endl;
	cout << "func compute: " << funcComp << endl;
}

//-----------------------------------------------------------------------------------------------
// ???????? ?????????
//-----------------------------------------------------------------------------------------------
int main()
{
	cout << scientific << setprecision(8);
	if (flag == 0)
		Fine();
	else
		Barrier();

	return 0;
}