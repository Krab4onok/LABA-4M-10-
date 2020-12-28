#include<iostream>
#include<cmath>
#include<locale.h>
#include<fstream>
#include<algorithm>
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
	setlocale(LC_ALL, "Russian");
	cout << "Нестационарное уравнение теплопроводности: " << endl;
	cout << "U't = 4U" << '"' << "+e^(-t)*sin(7*pi*x)+1" << " , X от 0 до 1, Т от 0 до 5" << endl;
	cout << "U(x,0) = 1 - x^2 - начальное условие(температура стержня в начальный момент времени)" << endl;
	cout << "U(0,t) = cos(t) - левое граничное условие 1-го рода(температура на левом торце стержня)" << endl;
	cout << "U(n,t) = sin(4*t) - правое граничное условие 1-го рода(температура на правом торце стержня)" << endl;
	cout << "Запишем неявную разностную схему: " << endl;
	cout << "(V(Xi,Tj)-V(Xi,Tj-1))/tou - 4 {V(Xi-1,Tj)-2*V(Xi,Tj)+V(Xi+1,Tj)}/h^2 = e^(-t)*sin(7*pi*Xi) +1" << endl;
	cout << "V(Xi,0) = 1 - (Xi)^2" << endl;
	cout << "V(0,Tj) = cos(Tj) " << endl;
	cout << "V(n,Tj) = sin(4*Tj) " << endl;
	cout << endl;
	cout << "V - численное решение, полученное с шагом h" << endl;
	cout << "V2 - численное решение, полученное с шагом h/2" << endl;
	cout << "max|v-v2| = E2" << endl;
	double X = 1.0;// ДИАПОЗОН ПО Х
	double T = 5.0;// ДИАПОЗОН ПО Т
	double gamma = 4.0;//gamma^2
	int n, m;// РАЗМЕРНОСТЬ РАЗБИЕНИЯ
	cout << "введите размерность сетки n= ";
	cin >> n;
	cout << "m= ";
	cin >> m;
	//n = 10;//РАЗБИЕНИЕ ПО Х
	//m = 10;//РАЗБИЕНИЕ ПО У
	double tou, h;
	tou = T / m;
	h = X / n;
	double eps = tou+(h*h);
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
			fi[i] = (exp(-1*j * tou) * sin(7 * M_PI * i * h) + 1)*tou + V[i][j-1];
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
	double MAX_E2 = -1000;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			MAX_E2 = max(MAX_E2, abs(error_rate(V[i][j], V2[2 * i][2 * j])));
		}
	}
	ofstream fout;
	fout.open("new.txt");
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
	int flag = 0;
		cout << "Точность численного решения(E2) = " << MAX_E2<<endl;
		cout << "Размерность сетки: " << "n= " << n << " m=" << m << endl;
		cout << "eps = tou+h^2= "<<eps<< endl;
		
		if (flag == 1)
	 { 
		 for (int i = -1; i < n+1; i++)
		 {
			if (i == -1)
			{
				cout << "№ Узла (i,j)" << '|';
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
				cout << endl;
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
					cout << j * tou << '|';
					cout.width(15);
					cout.setf(ios::left);
					cout << V[i][j] <</*scientific << */'|';
					cout.width(15);
					cout.setf(ios::left);
					cout << V2[2 * i][2 * j] << /*scientific << */'|';
					cout.width(15);
					cout.setf(ios::left);
					cout << abs(error_rate(V[i][j], V2[2 * i][2 * j])) <</*scientific << */'|';
					cout << endl;
				}
			}
	     }
		 cout << endl;
		 for (int j = 0; j < m + 1; j++)
		 {
			 cout.width(10);
			 cout.setf(ios::left);
			 cout << "Слой №" << j << '|';
			 for (int i = 0; i < n + 1; i++)
			 {
				 cout.width(10);
				 cout.setf(ios::left);
				 cout << V[i][j] << '|';
			 }
			 cout << endl;
		 }
	 } 
}