

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

float fs = 10000;//  b ma byc dzielnikiem fs
//float fs2 = 271;

void amplitude_modulation(vector<float> OX, vector<float> MT, vector<float> &ZA,float ka,float fn)
{
	float y;
	int i = 0;

	for (float x : MT)
	{
		y = (ka * MT[i] + 1) * cos(2 * PI_F * fn * OX[i]);
		ZA.push_back(y);
		i++;
	}


}

void phase_modulation(vector<float> OX, vector<float> MT, vector<float> &ZP, float kp,float fn)
{
	float y;
	int i = 0;

	for (float x : MT)
	{
		y = cos((2 * PI_F* fn * OX[i]) + (kp * MT[i]));
		ZP.push_back(y);
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

void function_p(vector<float> OX, vector<float> &OY, int N)
{
	int i = 0;
	float p = 0;
	float sum;

	for (float x : OX)
	{
		for (int n = 1; n <= N; n++)
		{
			sum = (cos(12 * OX[i] * pow(n, 2)) + cos(16 * OX[i] * n)) / pow(n, 2);

			p = p + sum;
		}
		OY.push_back(p);
		p = 0;
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
	/*
	zad3 

		1A Modulacja fazy: Fmin=998 Fmax=1002 W=4
		1A Modulacja amplitudy: Fmin=998 Fmax=1002 W=4

		1B Modulacja fazy: Fmin=994 Fmax=1007 W=13
		1B Modulacja amplitudy: Fmin=992 Fmax=1008 W=16

		1C Modulacja fazy: Fmin=935 Fmax=1085 W=150
		1C Modulacja amplitudy: Fmin=997 Fmax=1003 W=6

	*/

	vector<float> x1;
	vector<float> y1;
	vector<float> ZA1;
	vector<float> ZP1;

	vector<float> ZA2;
	vector<float> ZP2;

	vector<float> ZA3;
	vector<float> ZP3;

	
	DFT_Coeff DFTA;
	vector<float> ASAA;
	vector<float> FSAA;
	vector<float> AS_PAA;

	DFT_Coeff DFTAF;
	vector<float> ASAF;
	vector<float> FSAF;
	vector<float> AS_PAF;


	DFT_Coeff DFTB;
	vector<float> ASBA;
	vector<float> FSBA;
	vector<float> AS_PBA;

	DFT_Coeff DFTBF;
	vector<float> ASBF;
	vector<float> FSBF;
	vector<float> AS_PBF;


	DFT_Coeff DFTC;
	vector<float> ASCA;
	vector<float> FSCA;
	vector<float> AS_PCA;

	DFT_Coeff DFTCF;
	vector<float> ASCF;
	vector<float> FSCF;
	vector<float> AS_PCF;


	create_OX(0, 10000, 1 / fs, x1);
	//signal_tone(x1, y1);
	function_p(x1,y1,2);
	data_file2(x1, "x.txt");
	data_file2(y1, "y.txt");
	
	//--------zad1A--------
	amplitude_modulation(x1,y1,ZA1,0.5,1000);	
	data_file2(ZA1, "ZA1.txt");

	phase_modulation(x1, y1, ZP1, 0.5,1000);
	data_file2(ZP1, "ZP1.txt");
	//--------zad1B--------

	cout << x1.size()<<"\n";
	cout << y1.size() << "\n";
	cout << ZP1.size() << "\n";

	amplitude_modulation(x1, y1, ZA2, 6, 1000);
	data_file2(ZA2, "ZA2.txt");

	phase_modulation(x1, y1, ZP2, 2, 1000);
	data_file2(ZP2, "ZP2.txt");

	//--------zad1C--------
	amplitude_modulation(x1, y1, ZA3, 60, 1000);
	data_file2(ZA3, "ZA3.txt");

	phase_modulation(x1, y1, ZP3, 60, 1000);
	data_file2(ZP3, "ZP3.txt");

	//--------zad2--------

	DFT(ZA1, DFTA);
	amplitude_spectrum(DFTA, ASAA, y1);
	amplitude_spectrum_prime(ASAA, AS_PAA, y1);
	
	frequency_scale(DFTA, FSAA, fs);
	data_file2(FSAA, "FSAA.txt");
	data_file2(AS_PAA, "ASPAA.txt");




   
}


