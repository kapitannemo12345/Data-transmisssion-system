

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
float N = 1;
float fs = 1000;//TB*10 =fs !!!



typedef unsigned char BYTE;

string Binary_stream(int text)
{
	char buffer[33];
	_itoa_s(text, buffer, 2);
	string s(buffer);
	return s;
}

void information_signal(vector<float> &OY,string s)
{
	float y;
	int i = 0;
	float czas_trwania = Tb * fs;
	for (string::size_type i = 0; i < 10; i++) {//10 bo do 10 bitow ograniczyc
		
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

		if (IF[i] == 1)
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
	vector<float> x1;
	vector<float> y1;
	vector<float> IF1;

	int i = 'ab';
	string s;
	s=Binary_stream(i);
	cout << s << "\n";
	
	create_OX(0, 10000, 1 / fs, x1);
	information_signal(IF1, s);

	cout << "size OX:" << x1.size()<<"\n";
	cout << "size IF1:" << IF1.size()<<"\n";

	data_file2(x1, "x.txt");
	data_file2(IF1, "IF1.txt");

	amplitude_keying(x1, IF1, y1);
	data_file2(y1, "y1.txt");


	//int j = 'b';
	//char buffer[33]; //the variable you will store i's new value (binary value) in
	//char buffer2[33];
	
	//_itoa_s(i, buffer, 2);
	//_itoa_s(j, buffer2, 2);

	//printf("binary: %s\n", buffer);

	//string s(buffer);
	//string r(buffer2);

	//cout << s<<"\n";
	//cout << r << "\n";

	//string t ;
	//t += s;
	//t += r;
	//cout << t << "\n";
	
}

