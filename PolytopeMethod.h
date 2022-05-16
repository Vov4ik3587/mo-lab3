#pragma once
#include <stdio.h>
#include <vector>
#include <iostream>

using namespace std;

class PolytopeMethod
{
	const int N = 2,
		MAXITER = 100;
	const double EPS = 1E-7;
	vector<vector<double>> D; // ���������� ������ ���������
	vector<double> xMiddle, xNew;
	double fMin, fMax;
	int xMin, xMax;
	double alpha = 1,
		gamma = 2.5,
		beta = 0.5;

	
	double func(vector<double> x);

	int Prepare(vector<double>& x0);
	double FindMinMax();
	int FindMiddle();
	int Reflection();  // ���������
	int Expansion(double fTemp);   // ����������
	int Contraction(double fTemp); // ������
	int Shrink();      // ��������

public:
	int Algorithm(vector<double>& x0);
};