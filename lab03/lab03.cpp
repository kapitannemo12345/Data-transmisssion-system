

#include "pch.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

const float  PI_F = 3.14159265358979f;

float A = 1;
float B = 9;
float C = 4;

float fs = 10;//czestotliwosc 

void create_OX(float t_start, float t_end, float delta_t, vector<float> &v)
{
	for (float i = t_start; i <= t_end; i = i + delta_t)
	{
		v.push_back(i);
	}
}

class DFT_Coeff {
	public:
	//double real, img;
	vector<float> real;
	vector<float> img;
};

void DFT(vector<float> OY, DFT_Coeff &dft)
{
	int N = OY.size();
	int i = 0;
	float pr = 0;
	float pi = 0;
	float sum;

	for (float x : OY)
	{
		for (int n = 0; n < N; n++)
		{
			pr+= OY[n] * cos(-2*(PI_F* i * n));//k=i?
			pi+= OY[n] * sin(-2*(PI_F* i * n));			
		}
		pr = pr / N;
		pi = pi / N;

		dft.real.push_back(pr);
		dft.img.push_back(pi);

		pr = 0;
		pi= 0;
		i++;
	}

}

void amplitude_spectrum( DFT_Coeff dft, vector<float> &AS, vector<float> OY)//widmo amplitudowe
{
	int N = OY.size();
	float y;
	int i = 0;
	for (float x : dft.real)
	{
		y = sqrt( pow( dft.real[i] , 2 ) + pow(dft.img[i], 2));
		y = y * 2 / N;//skalowanie ?
		AS.push_back(y);
		i++;
	}
}

void frequency_scale(DFT_Coeff dft, vector<float> &FS)//skala częstotliwości
{
	float y;
	int i = 0;
	for (float x : dft.real)
	{
		y = i*( fs / dft.real.size() );// zmienna globalnaa=dft.real.size() d o poprawy
		FS.push_back(y);
		i++;
	}
}

void signal_tone(vector<float> OX, vector<float> &OY)
{
	//cout << "test";
	float y;
	int i = 0;

	for (float x : OX)
	{
		y = 1 * sin(2 * PI_F * B * OX[i] + C * PI_F);
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

int main()
{
	//----zad2-----

	vector<float> x1;
	vector<float> y1;
	DFT_Coeff DFT1;
	vector<float> AS1;
	vector<float> FS1;

	create_OX(0.0, 194, 0.1, x1);// nie za male fs? 0.1 0.001
	signal_tone(x1, y1);
	data_file(x1, y1, "function_s.txt");
	DFT(y1,DFT1);
	amplitude_spectrum(DFT1, AS1,y1);
	frequency_scale(DFT1, FS1);


	data_file(FS1, AS1, "aspec1.txt");
    

	DFT_Coeff dft_value;





}

