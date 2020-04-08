

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

void create_OX(float t_start, float t_end, float delta_t, vector<float> &v)
{
	//cout << "start"<<t_start;
	//cout << "end" << t_start;

	for (int  i = t_start; i <t_end; i++)
	{
		v.push_back((1/fs)*i);
	}
	//v.pop_back();
	

	cout << "test  ox:" << v.size()<<"\n";
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
			pr+= OY[n] * cos(2*(PI_F* i * n)/N);//k=i? -2=-k
			pi+= OY[n] * sin(-2*(PI_F* i * n)/N);			
		}
		//pr = pr / N;
		//pi = pi / N;

		dft.real.push_back(pr);
		dft.img.push_back(pi);

		pr = 0;
		pi= 0;
		i++;
	}

}


void I_DFT(DFT_Coeff dft, vector<float> &OY)
{
	int N = dft.real.size();
	int i = 0;
	float pr = 0;
	float pi = 0;
	float sum;

	for (float x : dft.real)
	{
		for (int n = 0; n < N; n++)
		{
			pr += ((dft.real[n] * cos(2 * (PI_F* i * n) / N) - dft.img[n] * sin(2 * (PI_F* i * n) / N))) / N;//k=i? -2=-k			;

		}

		//OY[i] = OY[i] / N;



		//pi = pi / N;

		OY.push_back(pr);


		pr = 0;

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

void frequency_scale(DFT_Coeff dft, vector<float> &FS,float fs)//skala częstotliwości
{
	float y;
	int i = 0;
	for (float x : dft.real)
	{
		y = i*( fs / dft.real.size() );// zmienna globalnaa=dft.real.size() d o poprawy fs =probkowanie
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

void quadratic_function(float a, float b, float c, vector<float> OX, vector<float> &OY)
{	
	float y;

	for (float x : OX)
	{
		y = a * x * x + b * x + c;
		OY.push_back(y);
	}
}

void function_y(vector<float> OX, vector<float> &OY, vector<float> QF)
{
	float y;
	int i = 0;

	for (float x : OX)
	{
		y = 2 * QF[i] * QF[i] + 12 * cos(x);
		OY.push_back(y);
		i++;
	}

}

void function_z(vector<float> OX, vector<float> &OY, vector<float> QF, vector<float> FY)
{
	float z;
	int i = 0;

	for (float x : OX)
	{
		z = sin(2 * PI_F * 7 * x)*QF[i] - 0.2 * log10(abs(FY[i]) + PI_F);
		OY.push_back(z);
		i++;
	}

}

void function_u(vector<float> OX, vector<float> &OY, vector<float> QF, vector<float> FY, vector<float> FZ)
{

	float u;
	int i = 0;

	for (float x : OX)
	{
		u = sqrt(abs(FY[i] * FY[i] * FZ[i])) - 1.8 * sin(0.4 * OX[i] * FZ[i] * QF[i]);
		OY.push_back(u);
		i++;
	}


}

void function_v(vector<float> OX, vector<float> &OY)
{

	float v;
	int i = 0;

	for (float x : OX)
	{
		if (0.22 > OX[i] && OX[i] >= 0.0)
		{
			v = (1 - 7 * OX[i]) * sin((2 * PI_F * OX[i] * 10) / (OX[i] + 0.04));

		}
		else if (0.22 <= OX[i] && OX[i] < 0.7)
		{
			v = 0.63 * OX[i] * sin(125 * OX[i]);

		}
		else if (1.0 >= OX[i] && OX[i] >= 0.7)
		{
			v = pow(OX[i], -0.662) + 0.77 * sin(8 * OX[i]);

		}

		OY.push_back(v);
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

int main()
{
	//----zad2-----

	vector<float> x1;
	vector<float> y1;	
	DFT_Coeff DFT1;	
	vector<float> AS1;	
	vector<float> FS1;
	vector<float> AS_P;
	
	vector<float> yi;
	DFT_Coeff I_DFT1;
	vector<float> AS_I;
	vector<float> FS_i;
	vector<float> AS_P_I;


	
	//-------zad2-------
	create_OX(0, 555, 1/fs, x1);//	
	signal_tone(x1, y1);
	data_file(x1, y1, "function_s.txt");
	data_file2(x1, "x.txt");
	data_file2(y1, "y.txt");
	DFT(y1,DFT1);
	amplitude_spectrum(DFT1, AS1,y1);
	amplitude_spectrum_prime(AS1, AS_P,y1);	
	frequency_scale(DFT1, FS1,fs);
	
	I_DFT(DFT1, yi);
	//amplitude_spectrum(I_DFT1, AS_I, y1);
	//amplitude_spectrum_prime(AS_I, AS_P_I, y1);
	//frequency_scale(I_DFT1, FS_i, fs);

	data_file2(x1, "xii.txt");
	data_file2(yi, "yii.txt");


	data_file2(FS1, "x.txt");
	data_file2(AS1, "y.txt");
	data_file2(FS1, "xp.txt");
	data_file2(AS_P, "yp.txt");
	//-------zad3-------
	vector<float> x2;
	vector<float> y_QF;
	DFT_Coeff DFT2;
	vector<float> AS2;
	vector<float> FS2;
	vector<float> AS_P2;
	create_OX(0, 555, 1 / fs, x2);
	quadratic_function(5, 5, 5, x2, y_QF);
	data_file2(x2, "xs.txt");
	data_file2(y_QF, "ys.txt");
	DFT(y_QF, DFT2);
	amplitude_spectrum(DFT2, AS2, y_QF);
	amplitude_spectrum_prime(AS2, AS_P2, y_QF);
	frequency_scale(DFT2, FS2,fs);
	data_file2(FS2, "sxp.txt");
	data_file2(AS_P2, "syp.txt");


	vector<float> y_FY;
	DFT_Coeff DFT3;
	vector<float> AS3;
	vector<float> FS3;
	vector<float> AS_P3;	
	function_y(x2, y_FY,y_QF);
	data_file2(x2, "xs.txt");
	data_file2(y_FY, "ys.txt");
	DFT(y_FY, DFT3);
	amplitude_spectrum(DFT3, AS3, y_QF);
	amplitude_spectrum_prime(AS3, AS_P3, y_QF);
	frequency_scale(DFT3, FS3, fs);
	data_file2(FS3, "sxp.txt");
	data_file2(AS_P3, "syp.txt");


	vector<float> y_FZ;
	DFT_Coeff DFT4;
	vector<float> AS4;
	vector<float> FS4;
	vector<float> AS_P4;
	function_z(x2, y_FZ, y_QF, y_FY);
	
	data_file2(x2, "xs.txt");
	data_file2(y_FZ, "ys.txt");
	DFT(y_FZ, DFT4);
	amplitude_spectrum(DFT4, AS4, y_QF);
	amplitude_spectrum_prime(AS4, AS_P4, y_QF);
	frequency_scale(DFT4, FS4, fs);
	data_file2(FS4, "sxp.txt");
	data_file2(AS_P4, "syp.txt");




	vector<float> y_FU;
	DFT_Coeff DFT5;
	vector<float> AS5;
	vector<float> FS5;
	vector<float> AS_P5;
	function_u(x2, y_FU, y_QF, y_FY, y_FZ);

	data_file2(x2, "xs.txt");
	data_file2(y_FU, "ys.txt");
	DFT(y_FU, DFT5);
	amplitude_spectrum(DFT5, AS5, y_QF);
	amplitude_spectrum_prime(AS5, AS_P5, y_QF);
	frequency_scale(DFT5, FS5, fs);
	data_file2(FS5, "sxp.txt");
	data_file2(AS_P5, "syp.txt");



	vector<float> y_FV;
	DFT_Coeff DFT6;
	vector<float> AS6;
	vector<float> FS6;
	vector<float> AS_P6;
	function_v(x2, y_FV);

	data_file2(x2, "xs.txt");
	data_file2(y_FV, "ys.txt");
	DFT(y_FV, DFT6);
	amplitude_spectrum(DFT6, AS6, y_QF);
	amplitude_spectrum_prime(AS6, AS_P6, y_QF);
	frequency_scale(DFT6, FS6, fs);
	data_file2(FS6, "sxp.txt");
	data_file2(AS_P6, "syp.txt");

	vector<float> y_FP2;
	DFT_Coeff DFT7;
	vector<float> AS7;
	vector<float> FS7;
	vector<float> AS_P7;
	function_p(x2, y_FP2, 55);

	data_file2(x2, "xs.txt");
	data_file2(y_FP2, "ys.txt");
	DFT(y_FP2, DFT7);
	amplitude_spectrum(DFT7, AS7, y_QF);
	amplitude_spectrum_prime(AS7, AS_P7, y_QF);
	frequency_scale(DFT7, FS7, fs);
	data_file2(FS7, "sxp.txt");
	data_file2(AS_P7, "syp.txt");


	
	//data_file(FS1, AS1, "aspec1.txt");
	//data_file(FS1, AS_P, "aspec2.txt");
    
	DFT_Coeff dft_value;

}

