

#include "pch.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

const float  PI_F = 3.14159265358979f;

float A = 5;
float B = 5;//4
float C = 5;

float fs = 185;//  b ma byc dzielnikiem fs
float fs2 = 271;

void amplitude_modulation(vector<float> OX, vector<float> MT, vector<float> &ZA,float ka)
{
	float y;
	int i = 0;

	for (float x : MT)
	{
		y = (ka*MT[i] + 1)*cos(2* PI_F*fs*OX[i]);
		ZA.push_back(y);
		i++;
	}


}

void phase_modulation(vector<float> OX, vector<float> MT, vector<float> &ZA, float kp)
{
	float y;
	int i = 0;

	for (float x : MT)
	{
		y = cos((2 * PI_F * OX[i]) + (kp * MT[i]));
		ZA.push_back(y);
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
			pr += OY[n] * cos(2 * (PI_F* i * n) / N);//k=i? -2=-k
			pi += OY[n] * sin(-2 * (PI_F* i * n) / N);
		}
		//pr = pr / N;
		//pi = pi / N;

		dft.real.push_back(pr);
		dft.img.push_back(pi);

		pr = 0;
		pi = 0;
		i++;
	}

}

void amplitude_spectrum(DFT_Coeff dft, vector<float> &AS, vector<float> OY)//widmo amplitudowe
{
	int N = OY.size();
	float y;
	int i = 0;
	for (float x : dft.real)
	{
		y = sqrt(pow(dft.real[i], 2) + pow(dft.img[i], 2));
		//y = y * 2 / N;//skalowanie ?
		AS.push_back(y);
		i++;
	}
}

void amplitude_spectrum_prime(vector<float> &AS, vector<float> &AS_P, vector<float> OY)//widmo amplitudowe
{


	int N = OY.size();
	float y;
	int i = 0;
	for (float x : AS)
	{


		y = 10 * log10(AS[i]);
		AS_P.push_back(y);


		//AS_P.push_back(y);
		i++;
	}
	/*
	for (int i = 0; i < AS.size(); i++)
	{
		if (AS_P[i] < threshold)
		{
			AS_P[i] = 0.0;
		}
	}
	*/
}

void frequency_scale(DFT_Coeff dft, vector<float> &FS, float fs)//skala częstotliwości
{
	float y;
	int i = 0;
	for (float x : dft.real)
	{
		y = i * (fs / dft.real.size());// zmienna globalnaa=dft.real.size() d o poprawy fs =probkowanie
		FS.push_back(y);
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
		myfile << OX[i] << endl;;
		i++;
	}

	myfile.close();

}

int main()
{
   
}


