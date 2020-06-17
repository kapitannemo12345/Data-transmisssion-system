#include "pch.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <bitset>
#include <ctime>
using namespace std;

int bit_amount =11;
int bit_amount_h = ((bit_amount * 8) / 4) * 7;

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
	//cout << "kod wjesciowy:" << "\n";

	for (int k = 0; k < X.size(); k++)//fs/4 czy inna???
	{
		//cout << X[k];

	}
	//cout << "\n";

	int d[4];
	int j = 0;

	for (int k = 0; k < (bit_amount * 8) / 4; k++)
	{

		for (int i = 0; i < 4; i++, j++)
		{
			d[i] = X[j];
		}

		//cout << "dane:" << "\n";


		for (int i = 0; i < 4; i++)
		{
			//cout << d[i];
		}

		//cout << "\n";
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
		//cout << "wektor KD= d*G:" << "\n";

		for (int i = 0; i < 7; i++)
		{
			//cout << KD[i];
		}
		//cout << "\n";


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

		//cout << "wektor S=KD*H: " << "\n";

		for (int i = 0; i < 3; i++)
		{
			//cout << S[i];
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

		//cout << "\n";
		//cout << "wektor KD:" << "\n";

		for (int i = 0; i < 7; i++)
		{
			//cout << KD[i];
		}

		//cout << "\n";
		//cout << "wektor KD_prime poprawiony:" << "\n";

		for (int i = 0; i < 7; i++)
		{
			cout << KD_p[i];
			Y.push_back(KD_p[i]);
		}
		//cout << "\n";
		//cout << "\n";

	}
}

void Hamming_decode(vector<int> &X, vector<int> &Y)
{
	const int R[4][7] = { {0,0,1,0,0,0,0},{0,0,0,0,1,0,0},{0,0,0,0,0,1,0},{0,0,0,0,0,0,1} };
	//cout << "HAMMING" << "\n";
	//cout << "kod wjesciowy do dekodowania:" << "\n";

	for (int k = 0; k < X.size(); k++)//fs/4 czy inna???
	{
		//cout << X[k];

	}
	//cout << "\n";

	int d[7];
	int j = 0;

	int limit = 22;

	for (int k = 0; k < (bit_amount * 8) / 4; k++)
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

		//cout << "dekodowane dane:" << "\n";

		for (int i = 0; i < 4; i++)
		{
			//cout << DEC[i];
			Y.push_back(DEC[i]);
		}
		//cout << "\n";


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


	czas_trwania = 1000;


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
	float f0 = (N + 1) / Tb;
	float f1 = (N + 2) / Tb;
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
	float f0 = (N + 1) / Tb;
	float f1 = (N + 2) / Tb;
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
	int period = 1000;//przedzial ilosc probek/ilosc bitow
	// przedzialy
	xp = 0;
	xk = period;
	n = 1000;//n BŁĄDD n =ilosc probek na tb czyli 100
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

	//cout << "size ox:" << OX.size();

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

	//cout << "size ox:" << OX.size();

	for (float x : OX)
	{
		myfile << OX[i] << endl;;
		i++;
	}

	myfile.close();

}

//SEKCJA4: dekodowanie sygnalu

void decode_TTL(vector<float> &TTL, vector<int> &s, int ograniczenie)//154 bity
{
	//cout << "rozmiar tego:"<< TTL.size()<<"\n";
	float y;
	int j = 0;
	float check = 1000;//to jest źle

	int change = 0;

	for (int i = 0; i < ograniczenie; i++)
	{
		for (j; j < check; j++)//zle?
		{
			if (TTL[j] == 1)
			{
				change = 1;

			}

		}

		if (change == 1)
		{
			change = 0;
			s.push_back(1);

		}
		else
		{
			s.push_back(0);
		}

		check = check + 1000;//1000 probek w ciagu 1tb
	}

}

void SB2S(vector<int> &X, string &s)
{
	string var = "";

	for (int i = 0; i < X.size(); i++)
	{
		if (X[i] == 0)
		{
			var += "0";
		}
		else
		{
			var += "1";
		}
	}

	//cout << var;

	string text = "";
	stringstream sstream(var);

	while (sstream.good())
	{
		std::bitset<8>bits;
		sstream >> bits;
		s += char(bits.to_ulong());
	}
}


//SEKCJA KONCOWA

void create_OX(float t_start, float t_end, float delta_t, vector<float> &v)
{
	//cout << "start"<<t_start;
	//cout << "end" << t_start;

	for (int i = t_start; i < t_end; i++)
	{
		v.push_back((1 / fs)*i);
	}
	//v.pop_back();


	//cout << "test  ox:" << v.size() << "\n";
}

//LAB11

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


void white_noise_add(vector<float> &F,float alfa)
{
	/* Generate a new random seed from system time - do this once in your constructor */
	srand(time(0));

	/* Setup constants */
	const static int q = 15;
	const static float c1 = (1 << q) - 1;
	const static float c2 = ((int)(c1 / 3)) + 1;
	const static float c3 = 1.f / c1;

	/* random number in range 0 - 1 not including 1 */
	float random = 0.f;

	/* the white noise */
	float noise = 0.f;

	int x;

	for (int i = 0; i < F.size(); i++)
	{
		random = ((float)rand() / (float)(RAND_MAX + 1));
		noise = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * c3;
		F[i] = F[i] * alfa + noise * (1 - alfa);
		
	}
}

void BER(vector<int> &F, vector<int> &F_noise)
{
	float total_bits = F.size();
	float error_bits = 0;
	for (int i = 0; i < F.size(); i++)
	{
		if (F[i] != F_noise[i])
		{
			error_bits++;
		}
	}

	float wynik = error_bits / total_bits;

	cout << "Ilosc blednych bitow: "<<error_bits << "\n";
	cout << "Wynik BER: " << wynik <<"\n";

}

int main()
{
	cout << "bazowy strumien binarny: " << "\n";
	string i = "ALA MA KOTA";//11*8=88 literek do hamminga petla wykonuje sie 88/4=22 //11*3=33
	string s;
	s = s;
	s = Binary_stream(i);
	cout << "hamming strumien binarny: " << "\n";
	
	//cout << s[14] << "\n";

	vector<int> dane;

	table(s, dane);
	//int k = 0;
	int range = 8;
	int  it = 0;
	dane.insert(dane.begin() + 24, it);//dodakowe 0 dla spacji
	dane.insert(dane.begin() + 48, it);//

	//cout << "rozmiaar: " <<X.size()<< "\n";

	for (int k = 0; k < dane.size(); k++)
	{
		//cout << dane[k];
	}

	//cout << "\n";

	vector<int> hamming;
	//vector<int> Z;

	Hamming(dane, hamming);

	//for (int k=0; k < hamming.size(); k++)
	//{
		//cout << hamming[k];
	//}
	//cout << "hamming przerobiony:" << "\n";
	//cout << "y:" << hamming.size()<< "\n";

	//Hamming(X, Y);

	//for (int k = 0; k < Y.size(); k++)
	//{
		//cout << Y[k];
	//}

	vector<float> IF1;//syngal informacyjny
	vector<float> x1;


	information_signal(IF1, hamming, bit_amount_h); // czemu 154 ad2
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
	
	
	//vector<float> ak1_N;
	//vector<float> fk1_N;
	//vector<float> pk1_N;

	white_noise_add(ak1,0.78);//0.78-(2 bledy) male ;0.74 sporadyczne (5 bledow);0.72-(12 bledy) czeste
	white_noise_add(fk1,0.16);//0.16-(1 blad) male;0.13 sporadyczne (5 bledow);0.10-(10 bledy) czeste
	white_noise_add(pk1,0.085);//0.085-(2 blady) male;sporadyczne 0.08-(6 bledow) ,0.05-(8 bledow) czeste

	//DFT_Coeff dft1;
	//vector<float> as_a;
	//vector<float> fs_a;
	//vector<float> asp_a;
	//DFT(ak1, dft1);
	//amplitude_spectrum(dft1, as_a, ak1);
	//amplitude_spectrum_prime(as_a, asp_a, ak1);
	//frequency_scale(dft1, fs_a, fs);
	//data_file2(fs_a, "fs_a.txt");//os ox
	//data_file2(asp_a, "asp_a.txt");//os ox

	//DFT_Coeff dft2;
	//vector<float> as_p;
	//vector<float> fs_p;
	//vector<float> asp_p;
	//DFT(pk1, dft2);
	//amplitude_spectrum(dft2, as_p, pk1);
	//amplitude_spectrum_prime(as_p, asp_p, pk1);
	//frequency_scale(dft2, fs_p, fs);
	//data_file2(fs_p, "fs_p.txt");//os ox
	//data_file2(asp_p, "asp_p.txt");//os ox


	//DFT_Coeff dft3;
	//vector<float> as_f;
	//vector<float> fs_f;
	//vector<float> asp_f;
	//DFT(fk1, dft3);
	//amplitude_spectrum(dft3, as_f, fk1);
	//amplitude_spectrum_prime(as_f, asp_f, fk1);
	//frequency_scale(dft3, fs_f, fs);
	//data_file2(fs_f, "fs_f.txt");//os ox
	//data_file2(asp_f, "asp_f.txt");//os ox



	
	X_TA(x1, ak1, AK_XA);

	vector<float> int1;
	integration(AK_XA, int1, bit_amount_h);// czestotliwosc czas trwania bitu
	data_file2(int1, "intA.txt");
	vector<float> MTA;
	treshold_mt(int1, MTA, 0.005);
	data_file2(MTA, "MTA.txt");

	vector<float> int3;
	vector<float> int4;
	vector<float> int_sum;
	vector<float> OY_F1_A;
	vector<float> OY_F1_B;
	X_TF(x1, fk1, OY_F1_A, OY_F1_B);


	integration(OY_F1_A, int3, bit_amount_h);// czestotliwosc / czas trwania bitu
	integration(OY_F1_B, int4, bit_amount_h);// czestotliwosc / czas trwania bitu
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
	data_file2(int_sum, "intF.txt");
	vector<float> MTF;
	treshold_mt(int_sum, MTF, 0.05);
	data_file2(MTF, "MTF.txt");

	vector<float> int2;
	vector<float> PK_XA;
	X_TP(x1, pk1, PK_XA);



	integration(PK_XA, int2, bit_amount_h);// 
	data_file2(int2, "intP.txt");
	vector<float> MTP;
	treshold_mt(int2, MTP, 0.005);
	data_file2(MTP, "MTP.txt");

	vector<int> d_ASK;
	vector<int> d_FSK;
	vector<int> d_PSK;

	decode_TTL(MTA, d_ASK, bit_amount_h);

	cout << "\n";

	cout << "kod hamminga dla ASK po demodulacji:" << "\n";
	for (int i = 0; i < d_ASK.size(); i++)
	{
		cout << d_ASK[i];
	}
	cout << "\n";
	decode_TTL(MTF, d_FSK, bit_amount_h);
	decode_TTL(MTP, d_PSK, bit_amount_h);
	cout << "kod hamminga dla FSK po demodulacji:" << "\n";
	for (int i = 0; i < d_FSK.size(); i++)
	{
		cout << d_FSK[i];
	}
	cout << "\n";

	cout << "kod hamminga dla PSK po demodulacji:" << "\n";
	for (int i = 0; i < d_PSK.size(); i++)
	{
		cout << d_PSK[i];
	}
	cout << "\n";
	vector<int> d_HA;
	vector<int> d_HF;
	vector<int> d_HP;


	Hamming_decode(d_ASK, d_HA);
	cout << "-------------------------------------" << "\n";
	Hamming_decode(d_FSK, d_HF);
	cout << "-------------------------------------" << "\n";
	Hamming_decode(d_PSK, d_HP);
	cout << "-------------------------------------" << "\n";


	cout << "kod po dekodowaniu kanalowym dla ASK:" << "\n";
	for (int i = 0; i < d_HA.size(); i++)
	{
		cout << d_HA[i];
	}
	cout << "\n";
	cout << "kod po dekodowaniu kanalowym dla FSK:" << "\n";
	for (int i = 0; i < d_HF.size(); i++)
	{
		cout << d_HF[i];
	}
	cout << "\n";
	cout << "kod po dekodowaniu kanalowym dla PSK:" << "\n";
	for (int i = 0; i < d_HP.size(); i++)
	{
		cout << d_HP[i];
	}

	string text1;
	string text2;
	string text3;

	SB2S(d_HA, text1);
	cout << "\n";

	SB2S(d_HF, text2);
	cout << "\n";

	SB2S(d_HP, text3);
	cout << "\n";
	cout << "napis po dekodowawniu strumienia binarnego dla ASK:" << "\n";
	cout << text1 << "\n";
	cout << "\n";
	cout << "napis po dekodowawniu strumienia binarnego dla PSK:" << "\n";
	cout << text3 << "\n";
	cout << "\n";
	cout << "napis po dekodowawniu strumienia binarnego dla FSK:" << "\n";
	cout << text2 << "\n";
	
	
	cout<< "\n";
	cout << "ASK:" << "\n";
	BER(dane,d_HA);
	cout << "FSK:" << "\n";
	BER(dane, d_HF);
	cout << "PSK:" << "\n";
	BER(dane, d_HP);
}