

#include "pch.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>

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
	for (auto i :text) {
		char bytes[8];
		_itoa_s(i, bytes, 2);
		output = output + "0" + string(bytes);//dodawanie 0 na koniec
	}
	cout << output << endl;
	return output;
}

void information_signal(vector<float> &OY,string s,int ograniczenie)
{
	float y;
	int i = 0;
	float czas_trwania = Tb * fs;//czas trwanie jednego bity to fs *tb
	for (string::size_type i = 0; i < ograniczenie; i++) {//10 bo do 10 bitow ograniczyc
		
		if (s[i]=='1')
		{
			//cout << s[i] << ' ';
			for (int j=0;j< czas_trwania;j++)
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
	float A1=1;
	float A2=8;

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

void signal_tone(vector<float> OX, vector<float> &OY,float A,float f,float FI)
{
	//cout << "test";
	float y;
	int i = 0;


	for (float x : OX)
	{
		
		y = A * sin(2 * PI_F * f * OX[i] + f * PI_F);
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

		kluczowanie amplitudy : Fmin = 1,75 Fmax = 2,25 W = 0.5
		
		kluczowanie czestotliwosci : Fmin = 2,875 Fmax = 3,25 W = 0,375
	
		kluczowanie fazy : Fmin = 1.688 Fmax = 2.375 W = 0,687
	
	*/

	//----zad3----
	vector<float> x1;
	vector<float> y1;
	vector<float> y2;
	vector<float> y3;

	vector<float> IF1;

	string i = "ab";
	string s;
	s =  s;
	s=Binary_stream(i);
	cout << s << "\n";
	
	

	create_OX(0, s.size()*fs, 1 / fs, x1);
	information_signal(IF1, s, s.size());

	cout << "size OX:" << x1.size()<<"\n";
	cout << "size IF1:" << IF1.size()<<"\n";

	data_file2(x1, "x.txt");
	data_file2(IF1, "IF1.txt");

	amplitude_keying(x1, IF1, y1);
	data_file2(y1, "y1.txt");

	frequency_keying(x1, IF1, y2);
	data_file2(y2, "y2.txt");

	phase_keying(x1, IF1, y3);
	data_file2(y3, "y3.txt");

	//----zad4----
	vector<float> x1_p;
	vector<float> y1_p;
	vector<float> y2_p;
	vector<float> y3_p;

	vector<float> IF_p1;

	create_OX(0, 1000, 1 / fs, x1_p);
	information_signal(IF_p1, s, 10);

	cout << "rozmiar size" << s.size()*fs << "\n";

	cout << "size OX_p:" << x1_p.size() << "\n";
	cout << "size IF_p1:" << IF_p1.size() << "\n";

	//data_file2(x1, "x.txt");
	//data_file2(IF1, "IF2.txt");

	amplitude_keying(x1_p, IF_p1, y1_p);
	//data_file2(y1, "y1.txt");

	frequency_keying(x1_p, IF_p1, y2_p);
	//data_file2(y2, "y2.txt");

	phase_keying(x1_p, IF_p1, y3_p);
	//data_file2(y3, "y3.txt");



	DFT_Coeff DFT1;
	vector<float> AS1;
	vector<float> FS1;
	vector<float> AS_P;
	
	DFT(y1_p, DFT1);
	amplitude_spectrum(DFT1, AS1, y1_p);
	amplitude_spectrum_prime(AS1, AS_P, y1_p);
	frequency_scale(DFT1, FS1, fs);
	data_file2(FS1, "xfs1.txt");
	data_file2(AS_P, "ASy.txt");


	DFT_Coeff DFT2;
	vector<float> AS2;
	vector<float> FS2;
	vector<float> AS_P2;

	DFT(y2_p, DFT2);
	amplitude_spectrum(DFT2, AS2, y2_p);
	amplitude_spectrum_prime(AS2, AS_P2, y2_p);
	frequency_scale(DFT2, FS2, fs);
	data_file2(FS2, "xfs2.txt");
	data_file2(AS_P2, "ASy2.txt");


	DFT_Coeff DFT3;
	vector<float> AS3;
	vector<float> FS3;
	vector<float> AS_P3;

	DFT(y3_p, DFT3);
	amplitude_spectrum(DFT3, AS3, y3_p);
	amplitude_spectrum_prime(AS3, AS_P3, y3_p);
	frequency_scale(DFT3, FS3, fs);
	data_file2(FS3, "xfs3.txt");
	data_file2(AS_P3, "ASy3.txt");
	
	
}

