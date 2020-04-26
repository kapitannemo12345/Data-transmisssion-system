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

void clock_singal(vector<float> OX, vector<float> &CLK, int limit)//fs wyznaczane prze ilos probek w ox???
{
	//int i = 0;

	float y = 1;//clk zaczyna sie od 1

	for (int j = 0; j < limit; j++)
	{
		for (int k = 0; k < fs / 2; k++)//fs/4 czy inna???
		{
			CLK.push_back(y);
		}

		if (y == 1)
		{
			y = 0;
		}
		else
		{
			y = 1;
		}

	}
}

void Manchester(vector<float> CLK,vector<float> TTL, vector<float> &MCH)//fs wyznaczane prze ilos probek w ox???
{

	float y=0;
	int i = 0;
	float fs_check = fs;

	for (int j = 0; j < fs/2; j++)//do 1 spadku (fs/2) clk manchester przyjmuje 0
	{
		cout << y << "\n";
		MCH.push_back(y);
		i++;

	}

	y = 1;//od 1 bo piersza pojdzie w dol po 2 ifie 1*-1

	for (i;i<bit_amount*fs;i++)
	{

		if (i ==fs_check)//sprawdzanie na wzroscie clk
		{
			if (TTL[i - 1] == TTL[i])
			{
				y = y * -1;

			}

			fs_check = fs_check + fs;
		}

		if (CLK[i - 1] >CLK[i])//spadek clk
		{
			y =y*-1;
		}

		MCH.push_back(y);
	}
		
}

void NRZI(vector<float> CLK, vector<float> TTL, vector<float> &NRZI)//fs wyznaczane prze ilos probek w ox???
{
	float y = 0;
	int i = 0;
	//float fs_check = fs;

	for (int j = 0; j < fs / 2; j++)//do 1 spadku (fs/2) clk manchester przyjmuje 0
	{
		cout << y << "\n";
		NRZI.push_back(y);
		i++;

	}

	y = -1;//od 1 bo piersza pojdzie w dol po 2 ifie 1*-1

	for (i; i < bit_amount*fs; i++)
	{
		if (CLK[i - 1] > CLK[i])//spadek clk
		{
			if (TTL[i] == 1)//zmiana stanu gdy ttl=1 PYTANIE?:czy NRZI zawsze powinine zejc 0 do -1 czy jest to uwarunkowane stanem TTS
			{
				y = y * -1;

			}
			else
			{

			}		
		}
		NRZI.push_back(y);
	}
	
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

int main()
{
	string i = "ab";
	string s;
	s = s;
	s = Binary_stream(i);
	cout << s << "\n";

	vector<float> x1;
	vector<float> TTL;
	vector<float> CLK;
	vector<float> MCH;
	vector<float> NRZI1;

	

	create_OX(0, s.size()*fs, 1 / fs, x1);
	information_signal(TTL, s, s.size());
	clock_singal(x1, CLK, 2 * s.size());
	Manchester(CLK, TTL, MCH);
	NRZI(CLK, TTL, NRZI1 );

	data_file2(x1, "x.txt");
	data_file2(TTL, "IF1.txt");
	data_file2(CLK, "CLK.txt");
	data_file2(MCH, "MCH.txt");
	data_file2(NRZI1, "NRZI.txt");
}