

#include "pch.h"
#include <iostream>

#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>

using namespace std;


//SEKCJA0: DANE

const float  PI_F = 3.14159265358979f;

float A = 5;
float B = 5;//4
float C = 5;

float Tb = 0.1;
float N = 2;
float fs = 1000;//TB*10 =fs !!!

typedef unsigned char BYTE;

//SEKCJA1: DANE pierwszy moduł

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

//SEKCJA2: KODOWANIE KANAŁOWE

void table(string text, vector<int> &X)
{
	int i = 0;

	for (i; i < text.size(); i++)
	{
		if (text[i] == '0')
		{
			X.push_back(0);
		}
		else
		{
			X.push_back(1);

		}
	}
}

void negate(int *X, int position)
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

	for (int k = 0; k < 22; k++)
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


		//KD[2] = 1; //test wykrywania błędu

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

		cout << "dekodowane dane:" << "\n";

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
		if (bit_error != 0)
		{
			//partity_bit=KD[7];
			negate(KD_p, bit_error - 1);
			int index;
			index = (KD_p[0] + KD_p[1] + KD_p[2] + KD_p[3] + KD_p[4] + KD_p[5] + KD_p[6]) % 2;//ponowne wyznaczenie bitu parzystasci po poprawie

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

		if (reject == 0)
		{
			cout << "\n";
			cout << "wektor KD_prime poprawiony:" << "\n";

			for (int i = 0; i < 8; i++)
			{
				cout << KD_p[i];
				Y.push_back(KD_p[i]);
			}
		}

		string status = "zachowany";

		if (reject == 1)
		{
			status = "odrzucony";

		}
		cout << "\n";
		cout << "status  pakietu:" << status << "\n";
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
			if (((j % 8) == 0) && (i = !0))//prezskkokk bitu parzystosci
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


//SEKCJA3A: AKS PSK FSK

void information_signal(vector<float> &OY, vector<int> &s, int ograniczenie)
{
	float y;
	int i = 0;
	float czas_trwania = Tb * fs;//czas trwanie jednego bity to fs *tb

	
	czas_trwania = 100;

	
	for (string::size_type i = 0; i < ograniczenie; i++) {//10 bo do 10 bitow ograniczyc


		//cout << "bit analizowany to:  " << s[i] << "\n";

		if (s[i] == 1)
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
	float A1 = 0;//wymagane 
	float A2 = 1;

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

//SEKCJA3b: DEMODULACJA AKS PSK FSK

void X_TA(vector<float> OX, vector<float> OY, vector<float> &OY_XA)
{
	int i = 0;
	float A2 = 1;
	float f1 = N / Tb;
	float y;

	for (float x : OX)
	{

		y = OY[i] * (A2 * sin(2 * PI_F * f1 * OX[i] + 1 * PI_F));
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
		y = OY[i] * (A1 * sin(2 * PI_F * f * OX[i] + fi1));
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

		y1 = OY[i] * (A1 * sin(2 * PI_F * f0 * OX[i] + 1 * PI_F));
		OY_XA1.push_back(y1);

		y2 = OY[i] * (A1 * sin(2 * PI_F * f1 * OX[i] + 1 * PI_F));
		OY_XA2.push_back(y2);


		y1 = 0;
		y2 = 0;
		i++;
	}
}

//SEKCJAkoncowa

void treshold_mt(vector<float> INT1, vector<float> &MT, float h)
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

void integration(vector<float> OX, vector<float> &INT1, int bit_amount)
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
	n = 100;//n BŁĄDD n =ilosc probek na tb czyli 100
	h = (xk - xp) / n;//h mozna wznazyc raz?
	//cout << "krok: h=" << h << endl;
	//cout << " h=" << INt.size() << endl;
	integration = 0;//calka
	int k = 0;
	int i = 0;
	int l = 0;
	for (int j = 0; j < bit_amount; j++)
	{

		for (k; k < n; k++)
		{
			integration += OX[k] * (1 / fs);//1 bit== fs dlaetgo dzialalo teraz jest problem
			INT1.push_back(integration);

		}

		/*
		for ( l; l < n; l++)
		{
			INT1.push_back(integration);
		}
		*/
		//cout << "Wynik calkowania metoda prostokatow: " << integration << endl;

		k = n;
		l = n;
		n = n + period;

		integration = 0;

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

void data_file3(vector<int> OX, string name)
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


int main()
{
	string i = "ALA MA KOTA";//11*8=88 literek do hamminga petla wykonuje sie 88/4=22
	string s;
	s = s;
	s = Binary_stream(i);
	cout << "bazowy strumien binarny: " << "\n";

	//cout << s[14] << "\n";

	
	vector<int> dane;

	table(s, dane);
	//int k = 0;
	int range = 8;
	int  it = 0;
	dane.insert(dane.begin() + 24, it);//dodakowe 0 dla spacji
	dane.insert(dane.begin() + 48, it);//

	//cout << "rozmiaar: " <<X.size()<< "\n";

	for (int k=0; k < dane.size(); k++)
	{
		cout << dane[k];
	}
			
	//cout << "\n";
	//cout << "\n";
	

	vector<int> hamming;
	//vector<int> Z;

	Hamming(dane,hamming);

	for (int k=0; k < hamming.size(); k++)
	{
		cout << hamming[k];
	}
	cout << "hamming przerobiony:" << "\n";
	cout << "y:" << hamming.size()<< "\n";

	//Hamming(X, Y);


	//for (int k = 0; k < Y.size(); k++)
	//{
		//cout << Y[k];
	//}
	
	vector<float> IF1;//syngal informacyjny
	vector<float> x1;

	
	information_signal(IF1, hamming, 154); // czemu 154 ad2
	create_OX(0, IF1.size(), 1 / fs, x1);//154 probki od IF1

	data_file2(x1, "x.txt");
	data_file2(IF1, "IF1.txt");

	//stad 3 sciezki psk fsk itd
	//MODULACJA

	vector<float> ak1;//amplitude keying	
	vector<float> fk1;//frequency keying
	vector<float> pk1;//phase keying

	amplitude_keying(x1, IF1, ak1);
	data_file2(ak1, "ak1.txt");

	frequency_keying(x1, IF1, fk1);
	data_file2(fk1, "fk1.txt");

	phase_keying(x1, IF1, pk1);
	data_file2(pk1, "pk1.txt");

	
	//DEMODULACJA
	vector<float> AK_XA;
	X_TA(x1, ak1, AK_XA);
	vector<float> int1;
	integration(AK_XA, int1, 154);// czestotliwosc / czas trwania bitu
	data_file2(int1, "intA.txt");
	vector<float> MTA;
	treshold_mt(int1, MTA, 0);
	data_file2(MTA, "MTA.txt");

	vector<float> int3;
	vector<float> int4;
	vector<float> OY_F1_A;
	vector<float> OY_F1_B;

	X_TF(x1,fk1 ,OY_F1_A, OY_F1_B);
	integration(OY_F1_A, int3, 154);// czestotliwosc / czas trwania bitu
	integration(OY_F1_B, int4, 154);// czestotliwosc / czas trwania bitu

	
	vector<float> int2;
	vector<float> PK_XA;
	X_TP(x1, pk1, PK_XA);	
	integration(PK_XA, int2, 154);// czestotliwosc / czas trwania bitu
	data_file2(int2, "intP.txt");
	vector<float> MTP;
	treshold_mt(int2, MTP, 0);
	data_file2(MTP, "MTP.txt");




}


