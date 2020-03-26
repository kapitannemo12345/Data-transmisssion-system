

#include "pch.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

const float  PI_F = 3.14159265358979f;

float A=1;
float B = 5;
float C = 5;

//10hz 1s/10= 0.1


void create_OX(float t_start, float t_end, float delta_t, vector<float> &v)
{
	for (float i = t_start; i <= t_end; i = i + delta_t)
	{
		v.push_back(i);
	}
}

void signal_tone(vector<float> OX, vector<float> &OY)
{
	cout << "test";
	float y;
	int i = 0;

	for (float x : OX)
	{
		y =A*sin(2 * PI_F * B * OX[i] + C*PI_F);
		OY.push_back(y);
		i++;
	}

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

	//cout << "size ox:" << OX.size();
	
	for (float x : OX)
	{
		myfile << OX[i] << ",";
		i++;
	}

	myfile.close();

}

void quantisation(vector<float> OY, vector<float> &Q,int q)
{
	int i = 0;
	//int y;
	float z;

	for (float y : OY)
	{
		z=floor(OY[i] *pow(2, q-1)+pow(2, q - 1)-0.5);// 2 ? w liczniku bo  2^16= 65 536 czyli wartosci  od 0 do 32 768 beda ujemne
		

		Q.push_back(z);
		i++;
	}

}

int main()
{

	//-----zad1-----
	vector<float> x1;
	vector<float> y1;
	create_OX(0.0, 5, 0.001, x1);// nie za male fs? 0.1 0.001
	signal_tone(x1, y1);
	data_file(x1, y1, "function_s.txt");
	data_file2(x1, "xs.txt");
	data_file2(y1, "ys.txt");
	//-----zad2-----
	vector<float> Q1;
	
	quantisation(y1, Q1,16);
	data_file(x1, Q1, "quantisation_s.txt");
	data_file2(Q1, "Qs_1.txt");
	data_file2(x1, "xs.txt");
	data_file2(Q1, "ys.txt");
	
	//-----zad3-----
	vector<float> Q2;
	vector<float> y2;
	vector<float> x2;

	create_OX(0.0, 5, 0.002, x2);// fs2 =0.5 fs1
	signal_tone(x2, y2);
	quantisation(y2, Q2, 8);
	data_file(x2, Q2, "quantisation_s_v2.txt");
	data_file2(Q1, "Qs_2.txt");
    
}

