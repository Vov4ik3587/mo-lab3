//#include <iostream>
//#include <fstream>
//#include <functional>
//#include <vector>
//#include <iomanip>
//
//using namespace std;
//
//vector<double> xk(2), dx(2);
//vector<double> S(2);
//double epsF = 1e-4;
//double epsN = epsF;
//const int maxiter = 1000;
//double f1, f2;
//double r; // Коэффициент штрафа
//int funcComp = 0; // Кол-во вычислений функции
//int flag = 1; // 1 для барьерных, 0 для штрафных
//
//// Заданная функция
//double func(double x, double y)
//{
//	return 4 * (y - x) * (y - x) + 3 * (x - 1) * (x - 1);
//}
//
//// Задача б (равенство)
//double func_h(double x, double y)
//{
//	return x + 1 - y;
//}
//
//// Задача а (неравенство)
//double func_g(double x, double y)
//{
//	return x + y + 1;
//}
//
//// Штрафные функции для задачи б
//double func_H(double x, double y)
//{
//	return abs(func_h(x, y));
//	//return func_h(x, y) * func_h(x, y);
//	//return pow(func_h(x, y), 10);
//}
//
//// Штрафные функции для задачи а
//double func_G(double x, double y)
//{
//	return 0.5 * (func_g(x, y) + abs(func_g(x, y)));
//	//return 0.5 * (func_g(x, y) + abs(func_g(x, y))) * 0.5 * (func_g(x, y) + abs(func_g(x, y)));
//	//return pow(0.5 * (func_g(x, y) + abs(func_g(x, y))), 10);
//}
//
//// Барьерные функции
//double func_Gb(double x, double y)
//{
//	//return -1. / func_g(x, y);
//	return -log(-func_g(x, y));
//}
//
//// Вспомогательная функция
//double GetQ(vector<double>& xk)
//{
//	// Задача а
//	return func(xk[0], xk[1]) + r * (/*func_H(xk[0], xk[1]) +*/ func_G(xk[0], xk[1]));
//	// Задача б
//	//return func(xk[0], xk[1]) + r * (func_H(xk[0], xk[1]) /*+ func_G(xk[0], xk[1])*/);
//}
//
//double GetFine(vector<double>& xk)
//{
//	// Задача а
//	return 1 * (/*func_H(xk[0], xk[1]) + */func_G(xk[0], xk[1])); //??????????
//	// Задача б
//	//return 1 * (func_H(xk[0], xk[1]) /*+ func_G(xk[0], xk[1])*/); //??????????
//}
//
//double GetQb(vector<double>& xk)
//{
//	return func(xk[0], xk[1]) + r * (func_Gb(xk[0], xk[1]));
//}
//
//double GetBarrier(vector<double>& xk)
//{
//	return 1 * (func_Gb(xk[0], xk[1])); // ?????????
//}
//
//vector<double> oper_minus(vector<double>& v1, vector<double>& v2)
//{
//	vector<double> res(v1.size());
//	for (int i = 0; i < v1.size(); i++)
//		res[i] = v1[i] - v2[i];
//
//	return res;
//}
//
//vector<double> oper_mult(vector<double>& v1, double c)
//{
//	vector<double> res(v1.size());
//	for (int i = 0; i < v1.size(); i++)
//		res[i] = v1[i] * c;
//
//	return res;
//}
//
//vector<double> oper_div(vector<double>& v1, double c)
//{
//	vector<double> res(v1.size());
//	for (int i = 0; i < v1.size(); i++)
//		res[i] = v1[i] / c;
//
//	return res;
//}
//
//void SearchIntervalWithMinimum(double& lyampred, double& lyamnext)
//{
//	vector<double> ltemp(2), ltemppred(2), ltempnext(2);
//	double lyam;
//	double deltaS = 0.001;
//	lyam = 1;
//
//	ltemp[0] = xk[0] + lyam * S[0];
//	ltemp[1] = xk[1] + lyam * S[1];
//
//	ltemppred[0] = xk[0] + (lyam - deltaS) * S[0];
//	ltemppred[1] = xk[1] + (lyam - deltaS) * S[1];
//
//	ltempnext[0] = xk[0] + (lyam + deltaS) * S[0];
//	ltempnext[1] = xk[1] + (lyam + deltaS) * S[1];
//
//	double f, fdp, fdm;
//	flag == 0 ? f = GetQ(ltemp) : f = GetQb(ltemp);
//	flag == 0 ? fdm = GetQ(ltemppred) : fdm = GetQb(ltemppred);
//	flag == 0 ? fdp = GetQ(ltempnext) : fdp = GetQb(ltempnext);
//	funcComp += 3;
//	double fnext;
//
//	if (f > fdp)
//	{
//		deltaS = deltaS;
//		fnext = fdp;
//	}
//
//	else if (f > fdm)
//	{
//		deltaS = -deltaS;
//		fnext = fdm;
//	}
//
//	else
//	{
//		lyampred = lyam - deltaS;
//		lyamnext = lyam + deltaS;
//		return;
//	}
//
//	lyamnext = lyam + deltaS;
//	double ft = f;
//
//	do
//	{
//		deltaS = 2 * deltaS;
//		lyampred = lyam;
//		lyam = lyamnext;
//		lyamnext = lyam + deltaS;
//
//		ltemp[0] = xk[0] + lyamnext * S[0];
//		ltemp[1] = xk[1] + lyamnext * S[1];
//
//		f = fnext;
//		flag == 0 ? fnext = GetQ(ltemp) : fnext = GetQb(ltemp);
//		funcComp++;
//
//	} while (fnext < f);
//
//	if (lyampred > lyamnext) swap(lyampred, lyamnext);
//}
//
//void GoldenRatio(double& a, double& b)
//{
//	vector<double>ltemp1(2), ltemp2(2);
//	double lyam1 = a + 0.381966011 * (b - a);
//	double lyam2 = b - 0.381966011 * (b - a);
//
//	ltemp1[0] = xk[0] + lyam1 * S[0];
//	ltemp1[1] = xk[1] + lyam1 * S[1];
//	double f1, f2;
//	flag == 0 ? f1 = GetQ(ltemp1) : f1 = GetQb(ltemp1);
//	funcComp++;
//	ltemp2[0] = xk[0] + lyam2 * S[0];
//	ltemp2[1] = xk[1] + lyam2 * S[1];
//	flag == 0 ? f2 = GetQ(ltemp2) : f2 = GetQb(ltemp2);
//	funcComp++;
//	int iter = 0;
//	while (abs(b - a) > epsN)
//	{
//		if (f1 < f2)//"root prinadlegit [a,lyam2]"
//		{
//			b = lyam2;
//			lyam2 = lyam1;
//			lyam1 = a + 0.381966011 * (b - a);
//			//x1 = a + b - x2;
//			f2 = f1;
//			ltemp1[0] = xk[0] + lyam1 * S[0];
//			ltemp1[1] = xk[1] + lyam1 * S[1];
//			flag == 0 ? f1 = GetQ(ltemp1) : f1 = GetQb(ltemp1);
//			funcComp++;
//		}
//		else //"root prinadlegit [lyam1,b]"
//		{
//			a = lyam1;
//			lyam1 = lyam2;
//			lyam2 = b - 0.381966011 * (b - a);
//			//x2 = a + b - x1;
//			f1 = f2;
//
//			ltemp2[0] = xk[0] + lyam2 * S[0];
//			ltemp2[1] = xk[1] + lyam2 * S[1];
//			flag == 0 ? f2 = GetQ(ltemp2) : f2 = GetQb(ltemp2);
//			funcComp++;
//		}
//	}
//}
//
//void Hook_G()// lyambda calculation using GoldenRatio method
//{
//	double fk;
//	flag == 0 ? fk = GetQ(xk) : fk = GetQb(xk);
//	funcComp++;
//	double epsH = epsF;
//
//	double fpred = 0;
//	double ftemp, ftemp2;
//	double df = abs(fk - fpred);
//	double a = 0.0, b = 0.0;
//	double lyam;
//	vector<double> xtemp(2);
//	vector<double> xknext(2);
//	int i;
//	for (i = 0; i < maxiter && df > epsH; i++)
//	{
//		dx[0] = 1;
//		dx[1] = 1;
//		int k;
//		dx = oper_mult(dx, 2);
//
//		fpred = fk;
//		do {
//			dx = oper_div(dx, 2);
//			k = 0;
//			xtemp = xk;
//
//			//исследующий поиск
//			//по х
//			xtemp[0] = xk[0] + dx[0];
//			flag == 0 ? ftemp = GetQ(xtemp) : ftemp = GetQb(xtemp);
//			funcComp++;
//			if (ftemp < fpred)
//			{
//				k++;
//			}
//			else
//			{ //больше или равно
//				xtemp[0] = xk[0] - dx[0];
//				flag == 0 ? ftemp = GetQ(xtemp) : ftemp = GetQb(xtemp);
//				funcComp++;
//				if (ftemp < fpred)
//				{
//					k++;
//				}
//				else
//				{
//					xtemp[0] = xk[0];
//					ftemp = fpred;
//				}
//			}
//			//по у
//			xtemp[1] = xk[1] + dx[1];
//			flag == 0 ? ftemp2 = GetQ(xtemp) : ftemp2 = GetQb(xtemp);
//			funcComp++;
//			if (ftemp2 < ftemp)
//			{
//				k++;
//			}
//			else
//			{
//				xtemp[1] = xk[1] - dx[1];
//				flag == 0 ? ftemp2 = GetQ(xtemp) : ftemp2 = GetQb(xtemp);
//				funcComp++;
//				if (ftemp2 < ftemp)
//				{
//					k++;
//				}
//				else
//				{
//					xtemp[1] = xk[1];
//					ftemp2 = ftemp;
//				}
//			}
//		} while (k == 0 && dx[0] > epsH&& dx[1] > epsH);
//
//		S = oper_minus(xtemp, xk);
//
//		SearchIntervalWithMinimum(a, b);
//		GoldenRatio(a, b);
//		double lyam = a;
//
//		xknext[0] = xk[0] + lyam * S[0];
//		xknext[1] = xk[1] + lyam * S[1];
//		if (flag == 1 && /*func_h(xknext[0], xknext[1]) != 0 &&*/ func_g(xknext[0], xknext[1]) > 0)
//		{
//		}
//		else
//			xk = xknext;
//		flag == 0 ? fk = GetQ(xk) : fk = GetQb(xk);
//		funcComp++;
//		df = abs(fk - fpred);
//	}
//}
//
//void Fine()
//{
//	double fine = GetFine(xk);
//
//	// Изменение стратегии
//	//r /= 2;
//	//r /= 100;
//	//r /= 50;
//	//r = sqrt(r);
//	//r = pow(r, 1 / 3);
//	//r *= r;
//	ofstream fout;
//	fout.open("out.txt");
//	fout.imbue(locale("Russian"));
//	fout << setprecision(15);
//	int ii;
//	cout << "fine=" << fine << " xk: " << xk[0] << " " << xk[1] << endl;
//	fout << "fine=" << fine << " xk: " << xk[0] << " " << xk[1] << endl;
//	for (ii = 0; ii < maxiter && abs(fine) > epsF; ii++)
//	{
//		//r *= 2;
//		//r *= 100;
//		//r *= 50;
//		r *= r;
//		//r = r * r*r;
//		Hook_G(); //найти новое значение xk
//		fine = GetFine(xk);
//		cout << "fine=" << fine << " xk: " << xk[0] << " " << xk[1] << endl;
//		fout << "fine=" << fine << " xk: " << xk[0] << " " << xk[1] << endl;
//	}
//
//	fout << ii << '\t' << funcComp << '\t' << func(xk[0], xk[1]) << '\t';
//	cout << "f = " << func(xk[0], xk[1]) << endl;
//	cout << "Fine Iter: " << ii << endl;
//	// Кол-во вычислений функции
//	cout << "func compute: " << funcComp << endl;
//}
//
//void Barrier()
//{
//	ofstream fout;
//	fout.open("out.txt");
//	fout.imbue(locale("Russian"));
//	fout << setprecision(15);
//	double barrier = 1;
//	double barrierpred = 0;
//	//r *= 2;
//	//r *= 1000;
//	//r =r*r*r;
//	int ii;
//	for (ii = 0; ii < maxiter && abs(barrier - barrierpred) > epsF; ii++)
//	{
//		barrierpred = barrier;
//		//r=pow(r,1./11);
//		r /= 12;
//		//r /= 1000;
//		Hook_G(); //найти новое значение xk
//
//		barrier = GetBarrier(xk);
//		cout << "bar=" << barrier << " xk: " << xk[0] << "\t" << xk[1] << "\t"<< r << "\t" << xk[0] + xk[1] << "\t" << func(xk[0],xk[1])<< endl;
//		fout << "bar=" << barrier << " xk: " << xk[0] << "\t" << xk[1] << endl;
//	}
//	fout << ii << '\t' << funcComp << '\t' << func(xk[0], xk[1]) << '\t';
//	cout << "f = " << func(xk[0], xk[1]) << endl;
//	cout << "Barrier Iter: " << ii << endl;
//	cout << "func compute: " << funcComp << endl;
//}
//
//int main()
//{
//	// Начальное приближение
//	xk[0] = -10;
//	xk[1] = -5;
//
//	// Начальная величина коэффициента штрафа
//	if (flag == 0)
//		r = 2;
//	else
//		r = 100;
//
//	cout << scientific << setprecision(6);
//	if (flag == 0)
//		Fine();
//	else
//		Barrier();
//
//	return 0;
//}