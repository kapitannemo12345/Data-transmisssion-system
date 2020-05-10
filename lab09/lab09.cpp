

#include "pch.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <sstream> 

using namespace std;

const float  PI_F = 3.14159265358979f;

float A = 5;
float B = 5;//4
float C = 5;

float Tb = 0.1;//1
float N = 2;//10
float fs = 10;//TB*10 =fs !!!
float bit_amount = 16;


typedef unsigned char BYTE;

string Binary_stream(string text)
{
	string output;
	for (auto i : text) {
		char bytes[8];
		_itoa_s(i, bytes, 2);
		output = output + "0" + string(bytes);//dodawanie 0 na koniec
	}
	//cout << output << endl;
	return output;
}

void table(string text,vector<int> &X)
{
	int i = 0;
	

	for ( i; i < text.size(); i++)
	{
		if (text[i]=='0')
		{
			X.push_back(0);
		}
		else
		{
			X.push_back(1);

		}
	}


}

void negate(int *X,int position)
{
	if (X[position] == 0)
	{

		X[position] = 1;
	}
	else
	{
		X[position] = 0;

	}
	
}

void Hamming(vector<int> &X, vector<int> &Y)
{
	cout << "kod wjesciowy:" << "\n";

	for (int k = 0; k < X.size(); k++)//fs/4 czy inna???
	{
		cout << X[k];
		
	}
	cout << "\n";

	int d[4];	
	int j = 0;
	
	for (int k = 0; k < 4; k++)
	{

		for (int i = 0; i < 4; i++, j++)
		{
			d[i] = X[j];
		}

		cout << "dane:" << "\n";
		

		for (int i = 0; i < 4; i++)
		{
			cout << d[i];
		}

		cout << "\n";
		//cout << "macierz G:" << "\n";

		const int GM[7][4] = { {1,1,0,1},{1,0,1,1},{1,0,0,0},{0,1,1,1},{0,1,0,0},{0,0,1,0},{0,0,0,1} };
		const int PM[3][7] = { {1,0,1,0,1,0,1},{0,1,1,0,0,1,1},{0,0,0,1,1,1,1} };
		/*
		for (int j = 0; j < 7; j++)
		{
			for (int i = 0; i < 4; i++)
			{
				cout << GM[j][i];

			}
			cout << "\n";
		}
		*/
		int KD[7];
		int sum = 0;

		for (int j = 0; j < 7; j++)
		{
			for (int i = 0; i < 4; i++)
			{
				sum = sum + GM[j][i] * d[i];
			}
			sum = sum % 2;
			KD[j] = sum;
			sum = 0;
		}

		//cout << "\n";
		cout << "wektor KD= d*G:" << "\n";

		for (int i = 0; i < 7; i++)
		{
			cout << KD[i];
		}
		cout << "\n";

		
		//cout << "macierz H:" << "\n";
		/*
		for (int j = 0; j < 3; j++)
		{
			for (int i = 0; i < 7; i++)
			{
				cout << PM[j][i];

			}
			cout << "\n";
		}
		*/

		int S[4];
		int KD_p[7];


		//KD[5] = 1; //test wykrywania błędu

		for (int i = 0; i < 7; i++)
		{
			KD_p[i] = KD[i];
		}

		
		sum = 0;

		for (int j = 0; j < 3; j++)
		{
			for (int i = 0; i < 7; i++)
			{
				sum = sum + (PM[j][i] * KD[i]);
				//cout << PM[j][i] << "\n";
				//cout << KD[i] << "\n";
				//cout << sum << "\n";
				//cout << "\n";
			}
			sum = sum % 2;
			S[j] = sum;
			//cout << "\n";
			//cout << sum << "\n";

			sum = 0;
		}

		cout << "wektor S=KD*H: " << "\n";

		for (int i = 0; i < 3; i++)
		{
			cout << S[i];
		}

		int bit_error = 1 * S[0] + 2 * S[1] + 4 * S[2];
		if (bit_error != 0)
		{
			negate(KD_p, bit_error - 1);
			/*
			if (KD_p[bit_error-1] == 0)//bit_error-1 bo indeksowaniwe od 0
			{
				KD_p[bit_error-1] = 1;
			}
			else
			{
				KD_p[bit_error-1] = 0;
			}
			*/
		}

		cout << "\n";
		cout << "wektor KD:" << "\n";

		for (int i = 0; i < 7; i++)
		{
			cout << KD[i];
		}

		cout << "\n";
		cout << "wektor KD_prime poprawiony:" << "\n";

		for (int i = 0; i < 7; i++)
		{
			cout << KD_p[i];
			Y.push_back(KD_p[i]);
		}
		cout << "\n";
		cout << "\n";

		

	}
}

//pm=h
//gm = g

int main()
{
	string i = "ab";
	string s;
	s = s;
	s = Binary_stream(i);
	cout << s << "\n";

	cout << s[14] << "\n";	
	vector<int> X;
	
	table(s, X);
	
	for (int k = 0; k < X.size(); k++)
	{
		cout << X[k] ;
	}
	cout << "\n";
	cout << "\n";
	cout << "\n";

	vector<int> Y;
	

	//Hamming(X);
	Hamming(X,Y);




}

