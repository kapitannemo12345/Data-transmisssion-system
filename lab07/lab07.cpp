

#include "pch.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>

using namespace std;

const float  PI_F = 3.14159265358979f;

float A = 5;
float B = 5;//4
float C = 5;

float Tb = 1;
float N = 2;
float fs = 100;//TB*10 =fs !!!

typedef unsigned char BYTE;

string Binary_stream(string text)
{

	string output;
	for (auto i : text) {
		char bytes[8];
		_itoa_s(i, bytes, 2);
		output = output + "0" + string(bytes);//dodawanie 0 na koniec
	}
	cout << output << endl;
	return output;
}

void information_signal(vector<float> &OY, string s, int ograniczenie)
{
	float y;
	int i = 0;
	float czas_trwania = Tb * fs;//czas trwanie jednego bity to fs *tb
	for (string::size_type i = 0; i < ograniczenie; i++) {//10 bo do 10 bitow ograniczyc

		if (s[i] == '1')
		{
			//cout << s[i] << ' ';
			for (int j = 0; j < czas_trwania; j++)
			{
				y = 1;
				OY.push_back(y);

			}
		}
		else
		{
			//cout << ' ';
			for (int j = 0; j < czas_trwania; j++)
			{
				y = 0;
				OY.push_back(y);

			}
		}
	}
}

void amplitude_keying(vector<float> OX, vector<float> IF, vector<float> &OY)
{

	float y;
	int i = 0;
	float f = N / Tb;
	float A1 = 1;
	float A2 = 8;

	for (float x : OX)
	{

		if (IF[i] == 0)
		{
			y = A1 * sin(2 * PI_F * f * OX[i] + 1 * PI_F);
			OY.push_back(y);
			i++;

		}
		else
		{
			y = A2 * sin(2 * PI_F * f * OX[i] + 1 * PI_F);
			OY.push_back(y);
			i++;

		}
	}

}

void frequency_keying(vector<float> OX, vector<float> IF, vector<float> &OY)
{
	float y;
	int i = 0;
	float f0 = N + 1 / Tb;
	float f1 = N + 2 / Tb;
	float A1 = 1;

	for (float x : OX)
	{

		if (IF[i] == 0)
		{
			y = A1 * sin(2 * PI_F * f0 * OX[i] + 1 * PI_F);
			OY.push_back(y);
			i++;

		}
		else
		{
			y = A1 * sin(2 * PI_F * f1 * OX[i] + 1 * PI_F);
			OY.push_back(y);
			i++;

		}
	}
}

void phase_keying(vector<float> OX, vector<float> IF, vector<float> &OY)
{

	float y;
	int i = 0;
	float f = N / Tb;
	float A1 = 1;
	float fi0 = 0;
	float fi1 = PI_F;

	for (float x : OX)
	{

		if (IF[i] == 0)
		{
			y = A1 * sin(2 * PI_F * f * OX[i] + fi0);
			OY.push_back(y);
			i++;

		}
		else
		{
			y = A1 * sin(2 * PI_F * f * OX[i] + fi1);
			OY.push_back(y);
			i++;

		}
	}

}

void create_OX(float t_start, float t_end, float delta_t, vector<float> &v)
{
	//cout << "start"<<t_start;
	//cout << "end" << t_start;

	for (int i = t_start; i < t_end; i++)
	{
		v.push_back((1 / fs)*i);
	}
	//v.pop_back();


	cout << "test  ox:" << v.size() << "\n";
}

void data_file(vector<float> OX, vector<float> OY, string name)
{
	ofstream myfile;
	myfile.open(name);
	int i = 0;

	cout << "size ox:" << OX.size();
	//cout << "size oy:" << OX.size()<<endl;

	for (float x : OX)
	{
		myfile << OX[i] << " " << OY[i] << endl;
		i++;
	}

	myfile.close();
}

void data_file2(vector<float> OX, string name)
{
	ofstream myfile;
	myfile.open(name);
	int i = 0;

	cout << "size ox:" << OX.size();

	for (float x : OX)
	{
		myfile << OX[i] << endl;;
		i++;
	}

	myfile.close();

}


void X_TA(vector<float> OX, vector<float> OY, vector<float> &OY_XA)
{
	int i = 0;
	float A2 = 8;
	float f1 = N / Tb;
	float y;

	for (float x : OX)
	{

		y= OY[i] * ( A2 * sin(2 * PI_F * f1 * OX[i] + 1 * PI_F) );
		OY_XA.push_back(y);
	
		i++;

	}

}

void X_TP(vector<float> OX, vector<float> OY, vector<float> &OY_XP)
{
	float y;
	int i = 0;
	float f = N / Tb;
	float A1 = 1;
	float fi1 = PI_F;

	for (float x : OX)
	{

		y = OY[i]*( A1 * sin(2 * PI_F * f * OX[i] + fi1));
		OY_XP.push_back(y);
		i++;
		
	}

}

void X_TF(vector<float> OX, vector<float> OY, vector<float> &OY_XA1, vector<float> &OY_XA2)
{
	float y1;
	float y2;
	int i = 0;
	float f0 = N + 1 / Tb;
	float f1 = N + 2 / Tb;
	float A1 = 1;

	for (float x : OX)
	{

		y1 = OY[i]*( A1 * sin(2 * PI_F * f0 * OX[i] + 1 * PI_F));
		OY_XA1.push_back(y1);
		
		y2 = OY[i]*( A1 * sin(2 * PI_F * f1 * OX[i] + 1 * PI_F));
		OY_XA2.push_back(y2);
	

		y1 = 0;
		y2 = 0;
		i++;
	}
}

void treshold_mt(vector<float> INT1, vector<float> &MT,float h)
{
	int i = 0;
	
	float y;

	for (float x : INT1)
	{

		if (INT1[i] > h)
		{
			MT.push_back(1);

		}
		else
		{
			MT.push_back(0);

		}
		

		i++;

	}


}


void integration(vector<float> OX, vector<float> &INT1,int bit_amount)
{
	float xp;
	float xk;
	float h;
	float integration;
	float n;
	int prob_amount = OX.size();//ilsoc probek
	int period = 100;//przedzial ilosc probek/ilosc bitow
	// przedzialy
	xp = 0;
	xk = period;	
	n = fs;//n powinno rownac sie czestotliwosci probkowania
	h = (xk - xp) / n;//h mozna wznazyc raz?
	cout << "krok: h=" << h << endl;
	//cout << " h=" << INt.size() << endl;
	integration = 0;//calka
	int k = 0;
	int i = 0;
	int l = 0;
	for (int j=0; j < bit_amount; j++)
	{

		for (k; k < n; k++)
		{
			integration += OX[k]*(1/fs);

		}

		for ( l; l < n; l++)
		{
			INT1.push_back(integration);
		}
		
		cout << "Wynik calkowania metoda prostokatow: " << integration << endl;
		
		k = n;
		l = n;
		n = n + period;

		integration = 0;
	
	}

}

int main()
{
	
	vector<float> x1;
	vector<float> y1;
	vector<float> y2;
	vector<float> y3;

	vector<float> IF1;
	
	string i = "ab";
	string s;
	s = s;
	s = Binary_stream(i);
	cout << s << "\n";


	create_OX(0, 1000, 1 / fs, x1);
	information_signal(IF1, s, 10);

	cout << "size OX:" << x1.size() << "\n";
	cout << "size IF1:" << IF1.size() << "\n";

	data_file2(x1, "x.txt");
	data_file2(IF1, "IF1.txt");

	amplitude_keying(x1, IF1, y1);
	data_file2(y1, "y1.txt");

	frequency_keying(x1, IF1, y2);
	data_file2(y2, "y2.txt");

	phase_keying(x1, IF1, y3);
	data_file2(y3, "y3.txt");	
	vector<float> int1;
	vector<float> OY_XA;
	X_TA(x1, y1, OY_XA);
	data_file2(OY_XA, "OY_XA.txt");
	integration(OY_XA, int1, 10);// czestotliwosc / czas trwania bitu
	data_file2(int1, "int1.txt");
	vector<float> MT1;
	treshold_mt(int1, MT1,20);
	data_file2(MT1, "MT1.txt");



	vector<float> int2;
	vector<float> OY_XP;
	X_TA(x1, y3, OY_XP);
	data_file2(OY_XP, "OY_XP.txt");
	integration(OY_XP, int2, 10);// czestotliwosc / czas trwania bitu
	data_file2(int2, "int2.txt");
	vector<float> MT2;
	treshold_mt(int2, MT2, 0);
	data_file2(MT2, "MT2.txt");








	vector<float> int3;
	vector<float> int4;
	vector<float> OY_F1;
	vector<float> OY_F2;

	X_TF(x1, y2, OY_F1, OY_F2);
	data_file2(OY_F1, "OY_XF1.txt");
	data_file2(OY_F2, "OY_XF2.txt");

	integration(OY_F1, int3, 10);// czestotliwosc / czas trwania bitu
	integration(OY_F2, int4, 10);// czestotliwosc / czas trwania bitu

	
	

	vector<float> int_sum;

	float y = 0;
	int j = 0;
	
	for (float x : int3)
	{
		y = int4[j] - int3[j];//p1-p0?
		//cout << int3[j] << int3[j] << "\n";
		int_sum.push_back(y);
		y = 0;
		j++;
	}

	data_file2(int_sum, "intS.txt");

	vector<float> MT3;
	treshold_mt(int_sum, MT3, 0);
	data_file2(MT3, "MT3.txt");

	//data_file2(int1, "IF1.txt");

}

