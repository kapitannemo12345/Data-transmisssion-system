

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


		KD[2] = 1; //test wykrywania błędu

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

void Hamming_decode(vector<int> &X, vector<int> &Y)
{
	const int R[4][7] = { {0,0,1,0,0,0,0},{0,0,0,0,1,0,0},{0,0,0,0,0,1,0},{0,0,0,0,0,0,1} };
	cout << "HAMMING" << "\n";
	cout <<  "kod wjesciowy do dekodowania:" << "\n";

	for (int k = 0; k < X.size(); k++)//fs/4 czy inna???
	{
		cout << X[k];

	}
	cout << "\n";

	int d[7];
	int j = 0;

	for (int k = 0; k < 4; k++)
	{

		for (int i = 0; i < 7; i++, j++)
		{
			d[i] = X[j];
		}
		
		int DEC[4];
		int sum = 0;

		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < 7; i++)
			{
				sum = sum + R[j][i] * d[i];
			}
			sum = sum % 2;
			DEC[j] = sum;
			sum = 0;
		}
		
		cout<<"dekodowane dane:" << "\n";

		for (int i = 0; i < 4; i++)
		{
			cout << DEC[i];
			Y.push_back(DEC[i]);
		}
		cout << "\n";


	}
}


void SECDED(vector<int> &X, vector<int> &Y)
{
	cout << "SECDED" << "\n";
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
		cout << "macierz G:" << "\n";

		const int GM[8][4] = { {1,1,0,1},{1,0,1,1},{1,0,0,0},{0,1,1,1},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,1,1,0} };
		const int PM[4][8] = { {1,0,1,0,1,0,1,0},{0,1,1,0,0,1,1,0},{0,0,0,1,1,1,1,0},{1,1,1,1,1,1,1,1} };
		
		for (int j = 0; j < 7; j++)
		{
			for (int i = 0; i < 4; i++)
			{
				cout << GM[j][i];

			}
			cout << "\n";
		}
		
		int KD[8];
		int sum = 0;

		for (int j = 0; j < 8; j++)
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

		for (int i = 0; i < 8; i++)
		{
			cout << KD[i];
		}
		cout << "\n";

		
		cout << "macierz H:" << "\n";
		
		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < 7; i++)
			{
				cout << PM[j][i];

			}
			cout << "\n";
		}
		

		int S[4];
		int KD_p[8];


		KD[2] = 1; //test wykrywania błędu
		KD[4] = 0;


		for (int i = 0; i < 8; i++)
		{
			KD_p[i] = KD[i];
		}

		sum = 0;

		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < 8; i++)
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

		int reject = 0;
		int partity_bit;

		int bit_error = 1 * S[0] + 2 * S[1] + 4 * S[2];
		if ( bit_error != 0 )
		{
			//partity_bit=KD[7];
			negate(KD_p, bit_error - 1);
			int index;
			index=(KD_p[0] + KD_p[1] + KD_p[2] + KD_p[3] + KD_p[4] + KD_p[5] + KD_p[6]) % 2;//ponowne wyznaczenie bitu parzystasci po poprawie

			//cout << "co jest kurwa" <<KD_p[7]<< KD_p[0] <<KD_p[1] << KD_p[2] << KD_p[3]<<  KD_p[4] << KD_p[5] << KD_p[6]<<"\n" ;
			if (KD_p[7] != index)
			{
				reject = 1;
			}
		}

		cout << "\n";
		cout << "wektor KD:" << "\n";

		for (int i = 0; i < 8; i++)
		{
			cout << KD[i];
		}

		cout << "\n";
		cout << "wektor KD_prime poprawiony:" << "\n";

		for (int i = 0; i < 8; i++)
		{
			cout << KD_p[i];
			Y.push_back(KD_p[i]);
		}

		string status="zachowany";

		if (reject == 1) 
		{
			status = "odrzucony";

		}
		cout << "\n";
		cout << "status  pakietu:" <<status<< "\n";
		cout << "\n";
		cout << "\n";

	}
}

void SECDED_decode(vector<int> &X, vector<int> &Y)
{
	const int R[4][7] = { {0,0,1,0,0,0,0},{0,0,0,0,1,0,0},{0,0,0,0,0,1,0},{0,0,0,0,0,0,1} };
	cout << "SECDED" << "\n";
	cout << "kod wjesciowy do dekodowania:" << "\n";

	for (int k = 0; k < X.size(); k++)//fs/4 czy inna???
	{
		cout << X[k];

	}
	cout << "\n";

	int d[7];
	int j = 0;

	for (int k = 0; k < 4; k++)
	{

		for (int i = 0; i < 7; i++, j++)
		{
			if (((j % 8) == 0)&&(i=!0))//prezskkokk bitu parzystosci
			{
				i--;

			}
			else
			{
				d[i] = X[j];//spagetti code  
				
			}
						
		}
		
		int DEC[4];
		int sum = 0;

		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < 7; i++)
			{
				sum = sum + R[j][i] * d[i];
			}
			sum = sum % 2;
			DEC[j] = sum;
			sum = 0;
		}

		cout << "dekodowane dane:" << "\n";

		for (int i = 0; i < 4; i++)
		{
			cout << DEC[i];
			Y.push_back(DEC[i]);
		}
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
	cout <<"bazowy strumien binarny: "<< s << "\n";

	//cout << s[14] << "\n";	
	vector<int> X;
	
	table(s, X);

	/*
	for (int k = 0; k < X.size(); k++)
	{
		cout << X[k] ;
	}
	*/
	cout << "\n";
	cout << "\n";
	cout << "\n";

	vector<int> Y;
	vector<int> Z;

	//Hamming(X);
	Hamming(X,Y);


	for (int k = 0; k < Y.size(); k++)
	{
		cout << Y[k];
	}
	cout << "\n";
	Hamming_decode(Y, Z);


	cout << "-----------------SECDED-----------------"<<"\n";

	vector<int> Y2;
	vector<int> Z2;
	vector<int> xp;
	//X[0] = 1;
	//X[1] = 0;
	//X[2] = 1;
	//X[3] = 0;
	//X[4] = 1;
	//X[5] = 1;
	//X[6] = 1;
	//X[7] = 1;



	SECDED(X, Y2);
	

	
	for (int k = 0; k < Y2.size(); k++)
	{
		cout << Y2[k];
	}
	cout << "\n";

	SECDED_decode(Y2, Z2);



}

