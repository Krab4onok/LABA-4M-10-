#include<iostream>
#include<cmath>
#include<locale.h>
#include<fstream>
#define M_PI       3.14159265358979323846
using namespace std;

///////////////////////////////////////////////////////////////////////////////////
/////////////////////ЧИСЛЕННОЕ РЕШЕНИЕ ПРОГОНКОЙ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
void solution(double* V, double gamma, int n, int m, double T, double X, double*fi)//метод прогонки(для jой строки)
{
	double tou, h;
	tou = T / m;
	h = X / n;
	double* A = new double[n + 1];
	double* C = new double[n + 1];
	double* B = new double[n + 1];
	double* alfa = new double[n + 1];
	double* beta = new double[n + 1];
	for (int i = 1; i < n ; i++)//j-й слой
	{
		A[i] = gamma * tou / (h * h);//поддиагональ
		B[i] = gamma * tou / (h * h);//наддиагональ
		C[i] = (1 + ((2.0*gamma * tou) / (h * h)));//главная диагональ
	}
	alfa[1] = 0;
	beta[1] = fi[0];
	for ( int i = 1; i < n; i++)//прямой ход прогонки
	{
		alfa[i + 1] = B[i] / (C[i] - A[i] * alfa[i]);//коэффициенты для нахождения решения
		beta[i + 1] = (fi[i] + A[i] * beta[i]) / (C[i] - A[i] * alfa[i]);//^
	}
	V[n] = fi[n];
	for (int i = n-1; i >= 0; i--)//обратный ход
	{
		V[i] = alfa[i + 1] * V[i + 1] + beta[i + 1];
	}
}

double error_rate(double v, double v2)
{
	return v - v2;
}


void main()
{
	double X = 1.0;// ДИАПОЗОН ПО Х
	double T = 5.0;// ДИАПОЗОН ПО Т
	double gamma = 4.0;//gamma^2
	int n, m;// РАЗМЕРНОСТЬ РАЗБИЕНИЯ
	n = 20;//РАЗБИЕНИЕ ПО Х
	m = 20;//РАЗБИЕНИЕ ПО У
	double tou, h;
	tou = T / m;
	h = X / n;
	double eps = tou;
	if (eps <  (h * h))
	{
		eps = (h * h);
	}
	double* fi = new double[n+1];//ПРАВАЯ ЧАСТЬ СЛАУ
	double** V = new double* [n+1];//ЧИСЛЕННОЕ РЕШЕНИЕ ЗАДАЧИ ПОЛУЧЕННОЕ ПРОГОНКОЙ
	for (int i = 0; i < n+1; i++)
	{
		V[i] = new double[m + 1];
	}//V[i][j]
	for (int i = 0; i < n + 1; i++)
	{
		V[i][0] = 1 - (i * h) * (i * h);
	}
	for (int j = 1; j < m+1; j++)
	{
		
		double* buf = new double[n + 1];
		fi[0] = cos(j * tou);
		fi[n] = sin(4 * j * tou);
		for (int i = 1; i < n; i++)
		{
			fi[i] = (exp(j * tou) * sin(7 * M_PI * i * h) + 1)*tou + V[i][j-1];
		}
		solution(buf, 4, n, m, T, X, fi);
		for (int k = 0; k < n + 1; k++)
		{
			V[k][j] = buf[k];
		}
		
	}
	
	/// <summary>///////////////////
	/// V2(Xi,Tj)///////////////////
	/// </summary>//////////////////
	int n2 = 2*n;
	int m2 = 2 * m;
	double h2 = X / n2;
	double tou2 = T / m2;
    double** V2 = new double* [n2 + 1];//ЧИСЛЕННОЕ РЕШЕНИЕ ЗАДАЧИ ПОЛУЧЕННОЕ ПРОГОНКОЙ ПРИ УСЛОВИИ n = 2n, m=2m
	double* fi2 = new double[n2 + 1];
	for (int i = 0; i < n2 + 1; i++)
	{
		V2[i] = new double[m2 + 1];
	}
	for (int i = 0; i < n2 + 1; i++)
	{
		V2[i][0] = 1 - (i * h2) * (i * h2);
	}
	for (int j = 1; j < m2 + 1; j++)
	{
		double* buf2 = new double[n2 + 1];
		fi2[0] = cos(j * tou2);
		fi2[n2] = sin(4 * j * tou2);
		for (int i = 1; i < n2; i++)
		{
			fi2[i] = (exp(j * tou2) * sin(7 * M_PI * i * h2) + 1)*tou2 + V2[i][j-1];
		}
		solution(buf2, 4, n2, m2, T, X, fi2);
		for (int k = 0; k < n2 + 1; k++)
		{
			V2[k][j] = buf2[k];
		}
		
	}
	ofstream fout;
	fout.open("file.txt");
	fout << "x, t, v, " <<n<<" , "<<m<< endl;
	for (int i = 0; i < n+1; i++)
	{
		double x = i * h;
		for (int j = 0; j < m + 1; j++)
		{
			double t = j * tou;
			fout << x << "," << t << "," << V[i][j] << endl;
		}
	}
	fout.close();
	/// <summary>
	/// ВЫВОД ТАБЛИЦЫ
	/// </summary>
	setlocale(LC_ALL, "Russian");
	for (int i = -1; i < (n + 1); i++)
	{
		if (i == -1)
		{
			cout << "№ Узла (i,j)"<< '|'; 
			cout.width(10);
			cout.setf(ios::left);
			cout << "Xi " << '|';
			cout.width(10);
			cout.setf(ios::left);
			cout << "Tj " << '|';
			cout.width(15);
			cout.setf(ios::left);
			cout << " V(Xi,Tj) " << '|';
			cout.width(15);
			cout.setf(ios::left);
			cout << "V2(Xi,Tj)" << '|';
			cout.width(15);
			cout.setf(ios::left);
			cout << "E2(V-V2)" << '|';
			cout.width(15);
			cout.setf(ios::left);
			cout << "Точность соблюдена?" << '|';
			cout.width(15);
			cout.setf(ios::left);
			cout<<"eps"<<'|'<< endl;
		}
		else
		{
			for (int j = 0; j < (m + 1); j++)
			{
				cout << "(" << i << "," << j << ")       " << '|';
				cout.width(10);
				cout.setf(ios::left);
				cout << i * h << '|';
				cout.width(10);
				cout.setf(ios::left);
				cout << j *tou << '|';
				cout.width(15);
				cout.setf(ios::left);
				cout << V[i][j]<</*scientific << */'|';
				cout.width(15);
				cout.setf(ios::left);
				cout<<V2[2*i][2*j] << /*scientific << */'|';
				cout.width(15);
				cout.setf(ios::left);
				cout << abs(error_rate(V[i][j], V2[2 * i][2*j]))<</*scientific << */'|';
				if (abs(error_rate(V[i][j], V2[2 * i][2*j])) <= eps)
				{
					cout.width(18);
					cout.setf(ios::left);
					cout << "да" << '|';
				}
				else
				{
					cout.width(18);
					cout.setf(ios::left);
					cout << "нет" << '|';
				}
				cout.width(15);
				cout.setf(ios::left);
				cout << eps << '|';
				cout << endl;
			}
		}
	}
}