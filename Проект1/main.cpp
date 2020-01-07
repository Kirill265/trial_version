#include <iostream> // ��� ����������� �����-������
#include <string>   // ��� ����������� ������ �� �������� 
#include <fstream>  // ��� ���������� � ������ ������ �� ����� � � ����
#include <iomanip>  // ��� setprecision()
#include <cmath>    // �������������� ���������

using namespace std;//����������� ������������ ����������� ����

double F[11] = { 0.0 }, Y0[30], Y1[30], X0[30], data1[14500][40], data2[14500][40], A[3][3], B[3][3], Matrix[3][3], t = 0, shag;
double aa[5], b[5], c[2], t_zad, X_t, t_border, tangaj_1st, ryskanie_1st, kren_1st, delta_teta, delta_psi;

const double pi = 3.141592;
unsigned int choise = 0;  // ���� ��� ����
unsigned int rocket = 0;  // ���� ��� ������ ������
unsigned int target = 0;  // ���� ��� ������ ����
unsigned int choise2 = 0; // ���� ��� ����
double F0[11] = { 0.0 };    // { ��, ��, �����, �����, del_v, del_n, del_el, ������ �����������, ����������, ���� ������� ����������, ���� ����}

double V0 = 350;					// �������� �� ����� �� �����
double tangaj0 = 45.0 * pi / 180.0; // ������
double ryskanie0 = 0.0;				// ��������
double kren0 = 0.0;					// ����
double mass0 = 980.0;				// �����
double alf_Tg[3] = { 0.0 }, bet_Tg[3] = { 0.0 }; // ���� ����� � �����
double n_Ya_toF, n_Ya_toS;			//���������� ��� ������ � ������ � ������ �����

								  //������� �� ����� ����� � �����������
double Ro = cos(ryskanie0 / 2)*cos(tangaj0 / 2)*cos(kren0 / 2) - sin(ryskanie0 / 2)*sin(tangaj0 / 2)*sin(kren0 / 2);
double Lambda = sin(ryskanie0 / 2)*sin(tangaj0 / 2)*cos(kren0 / 2) + cos(ryskanie0 / 2)*cos(tangaj0 / 2)*sin(kren0 / 2);
double Mu = sin(ryskanie0 / 2)*cos(tangaj0 / 2)*cos(kren0 / 2) + cos(ryskanie0 / 2)*sin(tangaj0 / 2)*sin(kren0 / 2);
double Nu = cos(ryskanie0 / 2)*sin(tangaj0 / 2)*cos(kren0 / 2) - sin(ryskanie0 / 2)*cos(tangaj0 / 2)*sin(kren0 / 2);

// ��� ����� � �����
//                  0             1                     2            3    4     5   6    7          8      9      10  11      12  13   14     15       16      17       18       19         20         21        22         23        24        25         26         27        28          29
//               {  ?,           Vxg,                  Vyg,         Vzg,  Xg,  Yg, Zg, Omega x, Omega y, Omega z, Ro, Lambda, Mu, Nu, kren, ryskanie, tangaj, �����, x ���� 1, y ���� 1, z ���� 1, Vx ���� 1, Vy ���� 1, Vz ���� 1, x ���� 2, y ���� 2, z ���� 2 , Vx ���� 2, Vy ���� 2, Vz ���� 2 }
double X00[30] = { 0.0, V0 * cos(tangaj0),	V0 * sin(tangaj0), 0.0, 0.0, 0.0, 0.0, 0.0,      0.0,     0.0,   Ro, Lambda, Mu, Nu, kren0, ryskanie0, tangaj0, mass0,  30000.0,	8000.0,	 0.0,     -350.0,		 0.0,    0.0,     33000.0, 7000.0,    0.0,     -400.0,       0.0,       0.0 };
 
//��� �������� ������ � 1-�� � 2-�� �����:
double Xf[30] = { 0.0 }, Xs[30] = { 0.0 }; //��������� ������ � ����������
double tf = 0, ts = 0; //����� ������
double F_M[11] = { 0.0 }; //���������� � ������ ��������� F
double dr_dt1[14500] = { 0.0 }, dr_dt2[14500] = { 0.0 }; //����������� �� ������� �������� ����������� �� ���� 1 � 2 ��������������
double var_H, var_L; //��������� ���������� ���� ���������

double dRo = 0.0;
double dLambda = 0.0;
double dMu = 0.0;
double dNu = 0.0;

double Tangaj = 0.0;
double Kren = 0.0;
double Ryskanie = 0.0;

//�������������� ���������� �������-����������
void IntRG(const double Ro, const double Lambda, const double Mu, const double Nu, const double OmegaX, const double OmegaY, const double OmegaZ)
{
	dRo = -(OmegaX*Lambda + OmegaY*Mu + OmegaZ*Nu) / 2;
	dLambda = (OmegaX*Ro - OmegaY*Nu + OmegaZ*Mu) / 2;
	dMu = (OmegaX*Nu + OmegaY*Ro - OmegaZ*Lambda) / 2;
	dNu = (-OmegaX*Mu + OmegaY*Lambda + OmegaZ*Ro) / 2;
}

//������ �� ������������ � ���� ������
void Ugl(const double Ro, const double Lambda, const double Mu, const double Nu)
{
	Ryskanie = atan2(2 * (Ro*Mu - Lambda*Nu), pow(Ro, 2) + pow(Lambda, 2) - pow(Mu, 2) - pow(Nu, 2));
	Kren = atan2(2 * (Ro*Lambda - Mu*Nu), pow(Ro, 2) - pow(Lambda, 2) + pow(Mu, 2) - pow(Nu, 2));
	Tangaj = asin(2 * (Ro*Nu + Mu*Lambda));
}

//����� ������������ �� �����
double min(double x1, double x2)
{
	if (x1 <= x2)
		return x1;
	else //x2 < x1
		return x2;
}
// ����� ������������ �� �����
double max(double x1, double x2)
{
	if (x1 >= x2)
		return x1;
	else //x2 > x1
		return x2;
}



// ��������� �������
double sk_zv(double height);        //������� ������� �������� ����� � ���������
double pl(double height);			//������� ������� ��������� ���������
double cx(double M, double alfa_);  //������� ���������� ��� ��������� (������� ����������� �� ����� �����, �����, ������)
double cy_a(double M, double alfa_);
double cz_b(double M, double alfa_);
double cy_d(double M, double alfa_);
double cz_d(double M, double alfa_);
double mx_wx(double M, double alfa_);
double mx_d(double M, double alfa_);
double mz_wz(double M, double alfa_);
double mz_a(double M, double alfa_);
double mz_d(double M, double alfa_);
double my_wy(double M, double alfa_);
double my_b(double M, double alfa_);
double my_d(double M, double alfa_);

double viz_fi(double Y, double r_viz, double Yc);        //������� ���������� ������� ����� ����������� � ������������ ���������
double viz_hi(double X, double Xc, double Z, double Zc); //������� ���������� ������� ����� ����������� � �������������� ���������

double r_vizV(double X, double Y, double Z, double Xc, double Yc, double Zc); //������� ���������� ������� ����� �����������

void Eyler(double tau, double X[30]);
void RG_sv_st(double Ro, double Lambda, double Mu, double Nu);
void A_sv_st(double t, double r, double k);
void B_sf_st(double fi, double hi);
void M_sv_sk(double alfa, double betta);

int main()
{
	setlocale(LC_ALL, "Russian");

	//���������� ��������� ������
	int sh = 0;
	int i_f = 0, i_s = 0;
	// ����� ���������� � ���� �������
	string X_names[38] = {
		"Vxg                ", "Vyg                ", "Vzg                ",
		"Xg                 ", "Yg                 ", "Zg                 ",
		"Omega x            ", "Omega y            ", "Omega z            ",
		"Ro                 ", "Lambda             ", "Mu                 ",
		"Nu                 ", "����               ", "��������           ",
		"������             ", "�����              ", "x ���� 1           ",
		"y ���� 1           ", "z ���� 1           ", "Vx ���� 1          ",
		"Vy ���� 1          ", "Vz ���� 1          ", "x ���� 2           ",
		"y ���� 2           ", "z ���� 2           ", "Vx ���� 2          ",
		"Vy ���� 2          ", "Vz ���� 2          ", "��                 ",
		"��                 ", "�����              ", "�����              ",
		"del_v              ", "del_n              ", "del_el             ",
		"������ ����������� ", "����������         "
	};
	shag = 0.01;

	do
	{
		t = 0;
		sh = 0;
		var_H = 0;
		var_L = 0;
		for (int i = 0; i < 30; i++)
		{
			X0[i] = X00[i];
			if (i >= 0 && i < 11)
			{
				F[i] = F0[i];
			}
		}
		do
		{   // ����� ������
			cout << "�������� ���� �� ���� �����. ������� 1 ��� 2.\n";
			cin >> rocket;
			if (rocket != 1 && rocket != 2)
				cout << "��� ���������� ������������ ���� ������!\n";
		} while (rocket != 1 && rocket != 2);

		if (rocket == 1)
		{
			tangaj_1st = tangaj0;
			ryskanie_1st = ryskanie0;
		}
		else // rocket == 2
		{
			do
			{   // ����� ����
				cout << "�������� ���� �� ���� �����. ������� 1 ��� 2.\n";
				cin >> target;
				if (target != 1 && target != 2 && target != 0)
					cout << "��� ���������� ������������ ���� ������!\n";
			} while (target != 1 && target != 2 && target != 0);

//			cout << "������� �������� �������: "; //���� ��������� ����� ������� 1-2 �����
//			cin >> t_zad;
			t_zad = 0;

			//������������ ��������� ��������� �����
			tangaj0 = 45.0 * pi / 180.0;
			X0[16] = tangaj0;						//��������� ���� ������� ��� ������ ������
			//���������� ��� ����� �����
			Ro = cos(ryskanie0 / 2)*cos(tangaj0 / 2)*cos(kren0 / 2) - sin(ryskanie0 / 2)*sin(tangaj0 / 2)*sin(kren0 / 2);
			Lambda = sin(ryskanie0 / 2)*sin(tangaj0 / 2)*cos(kren0 / 2) + cos(ryskanie0 / 2)*cos(tangaj0 / 2)*sin(kren0 / 2);
			Mu = sin(ryskanie0 / 2)*cos(tangaj0 / 2)*cos(kren0 / 2) + cos(ryskanie0 / 2)*sin(tangaj0 / 2)*sin(kren0 / 2);
			Nu = cos(ryskanie0 / 2)*sin(tangaj0 / 2)*cos(kren0 / 2) - sin(ryskanie0 / 2)*cos(tangaj0 / 2)*sin(kren0 / 2);
			X0[10] = Ro;
			X0[11] = Lambda;
			X0[12] = Mu;
			X0[13] = Nu;
			tangaj_1st = X0[16];
			ryskanie_1st = X0[15];
			kren_1st = X0[14];
			X0[18] += X0[21] * t_zad;
			X0[19] += X0[22] * t_zad;
			X0[20] += X0[23] * t_zad;
			X0[24] += X0[27] * t_zad;
			X0[25] += X0[28] * t_zad;
			X0[26] += X0[29] * t_zad;
		}

		do
		{
			if (t >= 3 && F[7] > 500 && rocket == 2)
			{
				if (n_Ya_toF < 15)
				{
					tf = t;
					for (int i = 0; i < 30; i++)
					{
						Xf[i] = X0[i];
						if (i >= 0 && i < 11)
						{
							F_M[i] = F[i];
						}
					}
					target = 1;
					i_f = -1;
					do
					{
						i_f += 1;
						Eyler(tf, Xf);
						for (int i = 1; i < 30; i++)
						{
							Xf[i] += shag*Y0[i];
						}
						tf += shag;
						if (n_Ya_toF > F[8])
							n_Ya_toF = F[8];
						if (i_f == 1)
							dr_dt1[sh] -= F[7];
						/*
						if (i_f == 500)
						{
							dr_dt1[sh] += F[7];
							dr_dt1[sh] /= ((i_f - 1) - 1);
						}
						*/
					} while ((F[7] >= 50) && (Xf[24] >= 1000) && (Xf[5] >= 0) && (Xf[4] <= Xf[18]) && (Xf[16] <= pi) && (tf <= 48));
					dr_dt1[sh] += F[7];
					dr_dt1[sh] /= (i_f - 1);

					for (int i = 0; i < 11; i++)
					{
					F[i] = F_M[i];
					}
				}

				if (n_Ya_toS < 15)
				{
					ts = t;
					for (int i = 0; i < 30; i++)
					{
						Xs[i] = X0[i];
						if (i >= 0 && i < 11)
						{
							F_M[i] = F[i];
						}
					}
					target = 2;
					i_s = -1;
					do
					{
						i_s += 1;
						Eyler(ts, Xs);
						for (int i = 1; i < 30; i++)
						{
							Xs[i] += shag*Y0[i];
						}
						ts += shag;
						if (n_Ya_toF > F[8])
							n_Ya_toS = F[8];
						if (i_s == 1)
							dr_dt2[sh] -= F[7];
						/*
						if (i_s == 500)
						{
							dr_dt2[sh] += F[7];
							dr_dt2[sh] /= ((i_s - 1) - 1);
						}
						*/
					} while ((F[7] >= 50) && (Xs[24] >= 1000) && (Xs[5] >= 0) && (Xs[4] <= Xs[24]) && (Xs[16] <= pi) && (ts <= 48));
					dr_dt2[sh] += F[7];
					dr_dt2[sh] /= (i_s - 1);

					for (int i = 0; i < 11; i++)
					{
					F[i] = F_M[i];
					}
				}
				if (dr_dt1[sh] == 0)
				{
					if (dr_dt2[sh] == 0) //��� � �����, ��� �� ����
					{
						cout << "\n\tAHTUNG!!!\n��������� ��������� ����� �� ���� ������������:\n� 1 ����:  " << n_Ya_toF << "\t�� 2 ����:  " << n_Ya_toS << "\n";
					}
					else //�� ������ ������
					{
						target = 2;
					}
				}
				else
				{
					if (dr_dt2[sh] == 0) //� ������ ������
					{
						target = 1;
					}
					else //������������ ��������� ���������
					{
						if (dr_dt2[sh] > dr_dt1[sh]) //��������� ���� ������� ������� ����������� (������� ����� ��������� ����� � ��� ��� ����)
						{
						//	var_H += ;
						//	var_L += ;
							target = 0;
						}
						else if (dr_dt2[sh] < dr_dt1[sh]) //��������� ���� ������� ������� ����������� (������� ����� ��������� ������ �� ��� ��� ����)
						{
						//	var_H -= ;
						//	var_L -= ;
							target = 0;
						}
						else //dr_dt2[sh] == dr_dt1[sh]   //��� �������!
						{
							target = 0;
						}
					}
				}
			}
			
			//��������������
			Eyler(t, X0); //��������� ������ ������
			for (int i = 1; i < 30; i++)
			{
				X0[i] += shag*Y0[i];
			}
//			cout << X0[4] << "\t"; cout << X0[5] << "\n";
			t = t + shag;
			//������ ������� ���������� � �������
			for (int i = 0; i < 30; i++)
			{
				data2[sh][i] = X0[i];
				if (i >= 0 && i<11)
					data2[sh][i + 30] = F[i];
			}
			if (rocket == 1)
			{
				X_t = X0[18]; //��������� ����
//				t_border = t; //����� ������� ������ ������
			}
			else if (target == 1) // rocket == 2
			{
				X_t = X0[18];
			}
			else // rocket == 2 && target == 2
			{
				X_t = X0[24];
			}
			if (F[7] > 10000)
			{
				t_border = t; //����� ������� ������ ������
			}
			sh = sh + 1;
			//		cout << X0[5] << "\n";

		} while (((F[7] >= 500 && rocket == 2) || (F[7] >= 50 && rocket == 1)) && X_t >= 1000 && X0[5] >= 0 && X0[4] <= X_t && X0[16] <= pi && t <= 148);


		cout << "�������� ��������� ������ � ����� : \n";
		//����� �������� ���������� ������ � �����
		for (int i = 1; i < 30; i++)
		{
			cout << X_names[i - 1] << " " << X0[i] << "\n";
		}
		cout << "r_viz               " << F[7] << "\n";
		cout << "t                   " << t << "\n";
		cout << "alf_Tg              " << alf_Tg[0] << "\n";
		cout << "bet_Tg              " << bet_Tg[0]  << "\n";
		cout << "bet_Tg              " << bet_Tg[0] << "\n";
		cout << "t_border            " << t_border << "\n";
		cout << "n_Ya1               " << n_Ya_toF << "\n";
		cout << "n_Ya2               " << n_Ya_toS << "\n";
		do
		{
			cout << "������� 1, ���� ������ ��������� ������, ����� - 0.\n";
			cin >> choise;
			if (choise != 0 && choise != 1)
				cout << "��� ���������� ������������ ���� ������!\n";
		} while (choise != 1 && choise != 0);

		if (choise == 0)
		{
			do
			{
				cout << "������ �� �������� ���������� �������� � ����, ���� ��, �� ������� 1, ����� - 0.\n";
				cin >> choise2;
				if (choise2 != 0 && choise2 != 1)
					cout << "��� ���������� ������������ ���� ������!\n";
			} while (choise2 != 1 && choise2 != 0);

			//������ ���������� � ����
			ofstream ffout("fine.txt");
			if (choise2 == 1)
			{
//				ffout << "t, c \t" << "V, �/� \t" << "tetta, ����\t" << "X, �\t" << "Y, �\t" << "Z, �\t" << "omegaZ, ���/c\t"
//					<< "m, kg\t" << "Xc_1, �\t" << "Yc_1, �\t" << "Zc_1, �\t" << "Xc_2, �\t" << "Yc_2, �\t" << "Zc_2, �\t" << "n_Ya, -\t"
//					<< "fi, ����\t" << "hi, ����\t" << "r_viz �\t" << "alfa, ����\t" << "betta, ����\t" << "delta_v, ����\t" << "delta_n, ����\n";
				for (int i = 0; i < sh; i += 10)
				{
					Ugl(data2[i][10], data2[i][11], data2[i][12], data2[i][13]);

					ffout /*<< fixed << setprecision(2) */ << i*shag << "\t";
					ffout /*<< fixed << setprecision(4)*/ << sqrt(pow(data2[i][1], 2) + pow(data2[i][2], 2) + pow(data2[i][3], 2)) << "\t";
					ffout << fixed << setprecision(4) << Tangaj * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][4] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][5] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][6] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][9] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][17] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][18] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][19] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][20] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][24] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][25] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][26] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][38] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][30] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][31] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][37] << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][32] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][33] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][34] * 180 / 3.141592 << "\t";
					ffout /*<< fixed << setprecision(4)*/ << data2[i][35] * 180 / 3.141592 << "\t";
					ffout << dr_dt1[i] << "\t";
					ffout << dr_dt2[i] << endl;

				}
				sh = sh - 1;
				Ugl(data2[sh][10], data2[sh][11], data2[sh][12], data2[sh][13]);
				ffout /*<< fixed << setprecision(2) */ << sh*shag << "\t";
				ffout /*<< fixed << setprecision(4)*/ << sqrt(pow(data2[sh][1], 2) + pow(data2[sh][2], 2) + pow(data2[sh][3], 2)) << "\t";
				ffout << fixed << setprecision(4) << Tangaj * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][4] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][5] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][6] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][9] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][17] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][18] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][19] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][20] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][24] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][25] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][26] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][38] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][30] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][31] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][37] << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][32] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][33] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][34] * 180 / 3.141592 << "\t";
				ffout /*<< fixed << setprecision(4)*/ << data2[sh][35] * 180 / 3.141592 << "\t";
				ffout << dr_dt1[sh] << "\t";
				ffout << dr_dt2[sh] << endl;
			}
			do
			{
				cout << "\n���������� ��� �����-�� �������?\n��� ������ �� ��������� ������� 0, ��� ���������� �������� ������� 1.\n";
				cin >> choise;
				if (choise != 0 && choise != 1)
					cout << "��� ���������� ������������ ���� ������!\n";
			} while (choise != 0 && choise != 1);
		}

	} while (choise == 1);

	return 0;
}
//������� ���������� ������ ������ ��
void Eyler(double tau, double X[30])
{
	double d = 0.5;
	double Ix = 640;
	double Iy = 2239;//2239; //640
	double Iz = 2239; //640
	double m = X[17];
	double g = 9.80665;
	double a = sk_zv(X[5]);
	double ro = pl(X[5]);
	double l = 6.3;
	double t_1st = 3.0;				//����������������� ������������ ������� ����������
	Y0[17] = -18.174 + 0.375*tau;		//����� ��������� ��������� �������
	if (tau >= 48)		//������� ����� ������� �������
		Y0[17] = 0;

	double P = 6230.0 * (-Y0[17]);

	double deltaAlf_Tg[3] = { 0, 0, 0 };
	double deltaBet_Tg[3] = { 0, 0, 0 };
	double Vx_newTg[3] = { 0, 0, 0 }, Vy_newTg[3] = { 0, 0, 0 }, Vz_newTg[3] = { 0, 0, 0 };

	double r_vizir1 = r_vizV(X[4], X[5], X[6], X[18], X[19], X[20]);
/*
	if ((rocket == 1) || (rocket == 2 && tau > t_border - t_zad)) //�������: ��� 1-�� ������ ��� ��� 2-�� (����� ������ ������ ������)
	{
		double teta_min[3] = { 0.00057, 0.00057, 0.0 };						//����������� ���������� ����� ����������� � ������������ ��������� ��� ������������ ����
		double psi_min[3] = { 9.0, 8.0, 0.0 };							//����������� ���������� ����� ����������� � �������������� ��������� ��� ������������ ����
		for (int i = 0; i < 3; i++)			//������� � �������
		{
			teta_min[i] *= (pi / 180.0);
			psi_min[i] *= (pi / 180.0);
		}
		double alf_max[3] = { 30.0, 30.0, 0.0 };							//����������� ��������� ������ ���������� ���� � ������������ ���������
		double alf_min[3] = { -30.0, -30.0, 0.0 };							//���������� ��������� ������ ���������� ���� � ������������ ���������
		double bet_max[3] = { 15.0, 15.0, 0.0 };							//����������� ��������� ������ ���������� ���� � �������������� ���������
		double bet_min[3] = { -15.0, -15.0, 0.0 };							//���������� ��������� ������ ���������� ���� � �������������� ���������
		double H_min_Tg[3] = { 3000.0, 3000.0, 0.0 };						//���������� ��������� ������ ������ ����
		double H_max_Tg[3] = { X00[19] + 1000.0, X00[25] + 1000.0, 0.0 };	//����������� ��������� ������ ������ ����
		double B_min_Tg[3] = { 800.0, 800.0, 0.0 };							//����������� ��������� ���������� ���� ����� 
		double B_max_Tg[3] = { -800.0, -800.0, 0.0 };						//����������� ��������� ���������� ���� ������ 

		double r_vizir = r_vizV(X[4], X[5], X[6], X[18], X[19], X[20]); //������ ����������� ������ ����
		if (r_vizir <= 12000.0) //������� ������ ������� ������ ����
		{
			//������ � ������������ ���������
			if (delta_teta > teta_min[0] && alf_Tg[0] < alf_max[0] - 1)
			{
				deltaAlf_Tg[0] = 0.01; //� ��������
			}
			else if (delta_teta < -teta_min[0] && alf_Tg[0] > alf_min[0] + 2)
			{
				deltaAlf_Tg[0] = -0.02; //� ��������
			}
			else
			{
				deltaAlf_Tg[0] = 0;	//� ��������
			}
			if (X[19] < H_min_Tg[0] && alf_Tg[0] < 0)
			{
				deltaAlf_Tg[0] = 0.002; //� ��������
			}
			else if (X[19] > H_max_Tg[0] && alf_Tg[0] > 0)
			{
				deltaAlf_Tg[0] = -0.006; //� ��������
			}
			deltaAlf_Tg[0] *= pi / 180;		//������� � �������
											//������ � �������������� ���������
			if (delta_psi > psi_min[0] && bet_Tg[0] > bet_min[0] + 2)
			{
				deltaBet_Tg[0] = -0.002; //� ��������
			}
			else if (delta_psi < -psi_min[0] && bet_Tg[0] < bet_max[0] - 2)
			{
				deltaBet_Tg[0] = 0.002; //� ��������
			}
			else
			{
				deltaBet_Tg[0] = 0;	//� ��������
			}

			if (X[20] < B_min_Tg[0] && bet_Tg[0] > 0)
			{
				deltaBet_Tg[0] = -0.006; //� ��������
			}
			else if (X[20] > B_max_Tg[0] && bet_Tg[0] < 0)
			{
				deltaBet_Tg[0] = 0.006; //� ��������
			}
			deltaBet_Tg[0] *= pi / 180;		//������� � �������
		}

*/
		//������� �������� �� ���������� �� � ���������
		M_sv_sk(deltaAlf_Tg[0], deltaBet_Tg[0]);
		//��������� ������� ���������� ����
		Vx_newTg[0] = Matrix[0][0] * X[21] + Matrix[0][1] * X[22] + Matrix[0][2] * X[23];
		Vy_newTg[0] = Matrix[1][0] * X[21] + Matrix[1][1] * X[22] + Matrix[1][2] * X[23];
		Vz_newTg[0] = Matrix[2][0] * X[21] + Matrix[2][1] * X[22] + Matrix[2][2] * X[23];
/*
		alf_Tg[0] += deltaAlf_Tg[0];
		bet_Tg[0] += deltaBet_Tg[0];

		r_vizir = r_vizV(X[4], X[5], X[6], X[24], X[25], X[26]);	 //������ ����������� ������ ����

		if (r_vizir <= 12000.0)	//������� ������ ������� ������ ����
		{
			//������ � ������������ ���������
			if (delta_teta > teta_min[1] && alf_Tg[1] < alf_max[1] - 2)
			{
				deltaAlf_Tg[1] = 0.002; //� ��������
			}
			else if (delta_teta < -teta_min[1] && alf_Tg[1] > alf_min[1] + 4)
			{
				deltaAlf_Tg[1] = -0.004; //� ��������
			}
			else
			{
				deltaAlf_Tg[1] = 0;	//� ��������
			}
			if (X[25] < H_min_Tg[1] && alf_Tg[1] < 0)
			{
				deltaAlf_Tg[1] = 0.004; //� ��������
			}
			else if (X[25] > H_max_Tg[1] && alf_Tg[1] > 0)
			{
				deltaAlf_Tg[1] = -0.012; //� ��������
			}
			deltaAlf_Tg[1] *= pi / 180;		//������� � �������

											//������ � �������������� ���������
			if (delta_psi > psi_min[1] && bet_Tg[1] > bet_min[1] + 4)
			{
				deltaBet_Tg[1] = -0.004; //� ��������
			}
			else if (delta_psi < -psi_min[1] && bet_Tg[1] < bet_max[1] - 4)
			{
				deltaBet_Tg[1] = 0.004; //� ��������
			}
			else
			{
				deltaBet_Tg[1] = 0;	//� ��������
			}

			if (X[26] < B_min_Tg[1] && bet_Tg[1] > 0)
			{
				deltaBet_Tg[1] = -0.008; //� ��������
			}
			else if (X[26] > B_max_Tg[1] && bet_Tg[1] < 0)
			{
				deltaBet_Tg[1] = 0.008; //� ��������
			}
			deltaBet_Tg[1] *= pi / 180;		//������� � �������
		}
*/
		//������� �������� �� ���������� �� � ���������
		M_sv_sk(deltaAlf_Tg[1], deltaBet_Tg[1]);
		//��������� ������� ���������� ����
		Vx_newTg[1] = Matrix[0][0] * X[27] + Matrix[0][1] * X[28] + Matrix[0][2] * X[29];
		Vy_newTg[1] = Matrix[1][0] * X[27] + Matrix[1][1] * X[28] + Matrix[1][2] * X[29];
		Vz_newTg[1] = Matrix[2][0] * X[27] + Matrix[2][1] * X[28] + Matrix[2][2] * X[29];
/*
		alf_Tg[1] += deltaAlf_Tg[1];
		bet_Tg[1] += deltaBet_Tg[1];

	}
	if ((rocket == 2) && (tau < t_border - t_zad)) //�������: ��� 2-�� ������ (����� ������ ������ ��� �� ����������)
	{
		//����� ������ � ����� �� ������ 
		int t_otn = int((tau + t_zad) / shag);
		for (int i = 18; i < 21; i++)
			X[i] = data2[t_otn - 1][i];
		for (int i = 21; i < 24; i++)
			X[i] = data2[t_otn][i];
		for (int i = 24; i < 27; i++)
			X[i] = data2[t_otn - 1][i];
		for (int i = 27; i < 30; i++)
			X[i] = data2[t_otn][i];
	}
*/

	Ugl(X[10], X[11], X[12], X[13]);

	//������������ ������������
	double L1 = sqrt(pow(X[10], 2) + pow(X[11], 2) + pow(X[12], 2) + pow(X[13], 2));
	for (int i = 10; i < 14; i++)
		X[i] /= L1;

	double tangaj = Tangaj;
	double ryskanie = Ryskanie;
	double kren = Kren;

	double Vxg = X[1];
	double Vyg = X[2];
	double Vzg = X[3];

	double V = sqrt(pow(Vxg, 2) + pow(Vyg, 2) + pow(Vzg, 2));
	double q = ro*pow(V, 2) / 2;
	double S = pi*pow(d, 2) / 4;
	double G = -m *g;
	double M = V / a;

	//������� ��������� �� ��������� �� � ���������
	RG_sv_st(X[10], X[11], X[12], X[13]);

	//������� ��������� ������ �� ��������� �� � ���������
	double Vx = A[0][0] * Vxg + A[0][1] * Vyg + A[0][2] * Vzg;
	double Vy = A[1][0] * Vxg + A[1][1] * Vyg + A[1][2] * Vzg;
	double Vz = A[2][0] * Vxg + A[2][1] * Vyg + A[2][2] * Vzg;

	double alfa = -atan2(Vy, Vx);
	double betta = asin(Vz / V);
	double Q = asin(Vyg / V);						//���� ������� ����������
	double Put = atan2(-Vzg, Vxg);					//���� ����
	double alfa1 = sqrt(pow(alfa, 2) + pow(betta, 2));

	//��� ��������� ����� ������� ����������
	delta_teta = Q - F[9];
	delta_psi = Put - F[10];
//	cout << delta_teta << "\n";
	F[2] = alfa;
	F[3] = betta;
	F[9] = Q;
	F[10] = Put;

	//������ ���
	double Cx = cx(M, alfa1);
	double Cy_a = cy_a(M, alfa);
	double Cz_b = cz_b(M, betta);

	double Cydelta = cy_d(M, alfa);
	double Czdelta = cz_d(M, betta);

	double mx_w = mx_wx(M, alfa1)*q*S*l;
	double my_w = my_wy(M, betta)*q*S*l;
	double mz_w = mz_wz(M, alfa)*q*S*l;

	double my_betta = my_b(M, betta)*q*S*l;
	double mz_alfa = mz_a(M, alfa)*q*S*l;

	double Mxdel = mx_d(M, alfa1)*q*S*l;
	double Mydel = my_d(M, betta)*q*S*l;
	double Mzdel = mz_d(M, alfa)*q*S*l;

	//��������� ������������ ���������� ��������� � ������������ � ������ ������
	
	aa[1] = -mz_w*l / (Iz*V);
	aa[2] = -mz_alfa / Iz;
	aa[3] = -Mzdel / Iz;
	aa[4] = (P + Cy_a*q*S) / (m*V);
	aa[5] = Cydelta*q*S / (m*V);

	b[1] = -my_w*l / (Iy*V);
	b[2] = -my_betta / Iy;
	b[3] = -Mydel / Iy;
	b[4] = (-P + Cz_b*q*S) / (m*V);
	b[5] = -Czdelta*q*S / (m*V);

	c[1] = -mx_w*l / (Ix*V);
	c[3] = -Mxdel / Ix;

	double ksiCC_v = 0.35;
	double KCC_v = 0.95;
	double ksiCC_n = 0.35;
	double KCC_n = 0.95;
	double ksiCC_el = 0.35;
	double TCC_el = 0.01;

	double K_v = (aa[2] * aa[5] - aa[3] * aa[4]) / (aa[2] + aa[1] * aa[4]);
	double T1_v = -aa[3] / (aa[3] * aa[4] - aa[2] * aa[5]);
	double T_v = 1 / sqrt(aa[2] + aa[1] * aa[4]);
	double ksi_v = (aa[1] + aa[4]) / (2 * sqrt(aa[2] + aa[1] * aa[4]));

	double K_n = (b[2] * b[5] - b[3] * b[4]) / (b[2] + b[1] * b[4]);
	double T1_n = -b[3] / (b[3] * b[4] - b[2] * b[5]);
	double T_n = 1 / sqrt(b[2] + b[1] * b[4]);
	double ksi_n = (b[1] + b[4]) / (2 * sqrt(b[2] + b[1] * b[4]));

	double K_el = c[3] / c[1];
	double T_el = 1 / c[1];

	double K2_v = -2 * T_v*(ksi_v*T1_v - pow(ksiCC_v, 2)*T_v - sqrt(pow(T_v, 2)*pow(ksiCC_v, 4) - 2 * ksi_v*pow(ksiCC_v, 2)*T_v*T1_v + pow(T1_v*ksiCC_v, 2))) / (K_v*pow(T1_v, 2));
	double K1_v = KCC_v*(1 + K_v*K2_v) / (K_v*pow(T1_v, 2));

	double K2_n = -2 * T_n*(ksi_n*T1_n - pow(ksiCC_n, 2)*T_n - sqrt(pow(T_n, 2)*pow(ksiCC_n, 4) - 2 * ksi_n*pow(ksiCC_n, 2)*T_n*T1_n + pow(T1_n*ksiCC_n, 2))) / (K_n*pow(T1_n, 2));
	double K1_n = KCC_n*(1 + K_n*K2_n) / (K_n*pow(T1_n, 2));

	double K1_el = (2 * ksiCC_el*T_el - TCC_el) / (K_el*TCC_el);
	double K2_el = T_el / (K_el*pow(TCC_el, 2));

	double K_fi = 200;
	double Kk_v = 1;

	double K_hi = 0;
	double Kk_n = 1;

	double Xc, Yc, Zc;
	double Vxc, Vyc, Vzc;

	//��������� ���� ��� ������ ������
	if (rocket == 1)
	{
		if (tau < t_1st) //�� ������� ���������
		{
			Xc = 30000*cos(tangaj_1st)*cos(ryskanie_1st);
			Yc = 30000*sin(tangaj_1st);
			Zc = 30000*cos(tangaj_1st)*sin(ryskanie_1st);
			Vxc = 0;
			Vyc = 0;
			Vzc = 0;
		}
		else if (r_vizir1 >= 6000) //�� ������ ������� ����������
		{
			Xc = X[18];
			Yc = X[19];
			Zc = X[20];
			Vxc = X[21];
			Vyc = X[22];
			Vzc = X[23];
		}
		else //�� ������� ��������
		{
			Xc = X[18];
			Yc = X[19];
			Zc = X[20];
			Vxc = X[21];
			Vyc = X[22];
			Vzc = X[23];
		}
	}
	// ������� ���������� ������� ����� ��������� �� ������ ���������� ����
	double Xc_12 = min(X[18], X[24]) - var_L;
	double Yc_12 = max(X[19], X[25]) + var_H;
	double Zc_12 = 0.5*(X[20] + X[26]);

	double Vxc_12 = 0.5*X[21] + 0.5*X[27];
	double Vyc_12 = 0.5*X[22] + 0.5*X[28];
	double Vzc_12 = 0.5*X[23] + 0.5*X[29];

	double R_viz_12 = r_vizV(X[4], X[5], X[6], Xc_12, Yc_12, Zc_12);
	//����� ����
	if (rocket == 2)
	{
		if (tau < t_1st) //�� ����������� ������� ����������
		{
			Xc = 30000*cos(tangaj_1st)*cos(ryskanie_1st);
			Yc = 30000*sin(tangaj_1st);
			Zc = 30000*cos(tangaj_1st)*sin(ryskanie_1st);
			Vxc = 0;
			Vyc = 0;
			Vzc = 0;
		}
		else if (target == 0)
//		else if (R_viz_12 >= 10000)
		//		else if (t <= t_border - t_zad) //�� ������ ������� ����������
		{
			Xc = Xc_12;
			Yc = Yc_12;
			Zc = Zc_12;
			Vxc = Vxc_12;
			Vyc = Vyc_12;
			Vzc = Vzc_12;
		}
		else //�� ������� ���������� � ���������� �����
			 //			if (R_viz_12 >= 30) 
		{
			if (target == 1)
			{
				Xc = X[18];
				Yc = X[19];
				Zc = X[20];
				Vxc = X[21];
				Vyc = X[22];
				Vzc = X[23];
			}
			if (target == 2)
			{
				Xc = X[24];
				Yc = X[25];
				Zc = X[26];
				Vxc = X[27];
				Vyc = X[28];
				Vzc = X[29];
			}
		}
	}
	F[7] = r_vizV(X[4], X[5], X[6], Xc, Yc, Zc);

	//���� ������� ���� �����������
	double fi_2 = viz_fi(X[5], F[7], Yc);
	double hi_2 = viz_hi(X[4], Xc, X[6], Zc);
	//������� �������� �� ��������� �� � �����������
	B_sf_st(fi_2, hi_2);
	//������� �������� ��������� ���� � ������ �� ��������� �� � �����������
	double Vr	= B[0][0] * (Vxc - Vxg) + B[0][1] * (Vyc - Vyg) + B[0][2] * (Vzc - Vzg);
	double Vfi	= B[1][0] * (Vxc - Vxg) + B[1][1] * (Vyc - Vyg) + B[1][2] * (Vzc - Vzg);
	double Vhi	= B[2][0] * (Vxc - Vxg) + B[2][1] * (Vyc - Vyg) + B[2][2] * (Vzc - Vzg);
	//����������� ����� ������� ����� �����������
	double d_fi = Vfi / F[7];
	double d_hi = -Vhi / (F[7] * cos(fi_2));
	F[0] = fi_2;
	F[1] = hi_2;

	//����������� ����� ������� ���������� ������
	double d_tet	= X[8] * sin(kren) + X[9] * cos(kren);
	double d_psi	= 1 / cos(tangaj)*(X[8] * cos(kren) - X[9] * sin(kren));
	double d_gamma	= X[8] - tan(tangaj)*(X[8] * cos(kren) - X[9] * sin(kren));

	//��������� ���������� ����� ������ ����������� � ������ ������ �������
	double del_v	= K1_v*K_fi*Kk_v*d_fi - K2_v*d_tet;
	double del_n	= K1_n*K_hi*Kk_n*d_hi - K2_n*d_psi;
	double del_el	= - K1_el*d_gamma - K2_el*kren;

	//������������� ����������� �� ���������� ����
	if (del_v > pi / 6)
		del_v = pi / 6;
	if (del_v < -pi / 6)
		del_v = -pi / 6;

	if (del_n > pi / 6)
		del_n = pi / 6;
	if (del_n < -pi / 6)
		del_n = -pi / 6;

	if (del_el > pi / 6)
		del_el = pi / 6;
	if (del_el < -pi / 6)
		del_el = -pi / 6;

	F[4] = del_v;
	F[5] = del_n;
	F[6] = del_el;

	//�� ����
	double FX = -Cx*q*S + P;
	double FY = Cy_a*q*S*alfa + Cydelta*q*S*del_v;
	double FZ = Cz_b*q*S*betta + Czdelta*q*S*del_n;
	//������� ��������� �� ��������� �� � ���������
	RG_sv_st(X[10], X[11], X[12], X[13]);
	//�������� �� ��� �� ��������� �� � ���������
	double Xg = A[0][0] * FX + A[1][0] * FY + A[2][0] * FZ;
	double Yg = A[0][1] * FX + A[1][1] * FY + A[2][1] * FZ;
	double Zg = A[0][2] * FX + A[1][2] * FY + A[2][2] * FZ;

	double Ya = FX*sin(alfa) + FY*cos(alfa);
	double n_Ya = (P * sin(alfa) + Ya) / (m*g);	//���������� ����������
	F[8] = n_Ya;
	//�� �������
	double mx = mx_w*l*X[7] / V + Mxdel*del_el;
	double my = my_betta*betta + my_w*l*X[8] / V + Mydel*del_n;
	double mz = mz_alfa*alfa + mz_w*l*X[9] / V + Mzdel*del_v;

	//������ ������ ������ ���������� �������-����������
	IntRG(X[10], X[11], X[12], X[13], X[7], X[8], X[9]);
	//������ ������ ������ ��
	Y0[1] = Xg / m; 		                                //Vxg
	Y0[2] = (Yg + G) / m;	                             	//Vyg
	Y0[3] = Zg / m;		                            		//Vzg
	Y0[4] = X[1];		                              		//Xg
	Y0[5] = X[2];		                              		//Yg
	Y0[6] = X[3];		                                	//Zg
	Y0[7] = mx / Ix + (Iy - Iz)*X[8] * X[9] / Ix;					//Omega x
	Y0[8] = my / Iy + (Iz - Ix)*X[7] * X[9] / Iy;					//Omega y
	Y0[9] = mz / Iz + (Ix - Iy)*X[7] * X[8] / Iz;					//Omega z

	Y0[10] = dRo;
	Y0[11] = dLambda;
	Y0[12] = dMu;
	Y0[13] = dNu;

	Y0[18] = Vx_newTg[0];
	Y0[19] = Vy_newTg[0];
	Y0[20] = Vz_newTg[0];
	Y0[24] = Vx_newTg[1];
	Y0[25] = Vy_newTg[1];
	Y0[26] = Vz_newTg[1];

	Y0[21] = 0;
	Y0[22] = 0;
	
//	cout << Xc << "\t"; cout << Yc << "\t"; cout << t << "\t"; cout << M << "\n";
	
	/*
	Y0[14] = d_gamma;
	Y0[15] = d_psi;
	Y0[16] = d_tet;
	*/
}
//������� ��������� �� ��������� �� � ��������� ����� ��������� �-�
void RG_sv_st(double Ro, double Lambda, double Mu, double Nu)
{
	A[0][0] = Ro*Ro + Lambda*Lambda - Mu*Mu - Nu*Nu;
	A[0][1] = 2 * (Ro*Nu + Lambda*Mu);
	A[0][2] = 2 * (-Ro*Mu + Lambda*Nu);

	A[1][0] = 2 * (-Ro*Nu + Lambda*Mu);
	A[1][1] = Ro*Ro + Mu*Mu - Nu*Nu - Lambda*Lambda;
	A[1][2] = 2 * (Ro*Lambda + Nu*Mu);

	A[2][0] = 2 * (Ro*Mu + Nu*Lambda);
	A[2][1] = 2 * (-Ro*Lambda + Nu*Mu);
	A[2][2] = Ro*Ro + Nu*Nu - Lambda*Lambda - Mu*Mu;
}
//������� ��������� �� ��������� �� � ��������� ����� ���� ������
void A_sv_st(double t, double r, double k)
{
	A[0][0] = cos(t)*cos(r);
	A[0][1] = sin(t);
	A[0][2] = -cos(t)*sin(r);

	A[1][0] = -sin(t)*cos(r)*cos(k) + sin(r)*sin(k);
	A[1][1] = cos(t)*cos(k);
	A[1][2] = sin(t)*sin(r)*cos(k) + cos(r)*sin(k);

	A[2][0] = sin(t)*cos(r)*sin(k) + sin(r)*cos(k);
	A[2][1] = -cos(t)*sin(k);
	A[2][2] = -sin(t)*sin(r)*sin(k) + cos(r)*cos(k);
}
//������� �������� �� ��������� �� � �����������
void B_sf_st(double fi, double hi)
{
	B[0][0] = cos(fi)*cos(hi);
	B[0][1] = sin(fi);
	B[0][2] = -cos(fi)*sin(hi);
	B[1][0] = -sin(fi)*cos(hi);
	B[1][1] = cos(fi);
	B[1][2] = sin(fi)*sin(hi);
	B[2][0] = sin(hi);
	B[2][1] = 0;
	B[2][2] = cos(hi);
}
//������� �������� �� ���������� �� � ���������
void M_sv_sk(double alfa, double betta)
{
	Matrix[0][0] = cos(alfa)*cos(betta);
	Matrix[0][1] = -sin(alfa)*cos(betta);
	Matrix[0][2] = sin(betta);
	Matrix[1][0] = sin(alfa);
	Matrix[1][1] = cos(alfa);
	Matrix[1][2] = 0;
	Matrix[2][0] = -cos(alfa)*sin(betta);
	Matrix[2][1] = sin(alfa)*sin(betta);
	Matrix[2][2] = cos(betta);
};
//������� ������� �������� ����� � ���������
double sk_zv(double H)
{
	double SoundSpeed;
	if (H <= 0.0)
	{
		SoundSpeed = 340.0;
	}
	else if (H > 0.0 && H <= 11000.0)
	{
		SoundSpeed = sqrt(299.79e-9*H*H - 2.6135*H + 115.80e+3);
	}
	else if (H > 11000.0 && H <= 25000.0)
	{
		SoundSpeed = 295.07;
	}
	else if (H > 25000.0 && H <= 46000.0)
	{
		SoundSpeed = sqrt(-179.80e-9*H*H + 1.1098*H + 59.433e+3);
	}
	else
	{
		SoundSpeed = 331.82;
	}
	return SoundSpeed;
};
//������� ������� ��������� ���������
double pl(double H)
{
	double DensutyAir;
	if (H <= 0.0)
	{
		DensutyAir = 1.225;
	}
	else if (H > 0.0 && H <= 13000.0)
	{
		DensutyAir = 3.1088e-9*H*H - 114.04e-6*H + 1.223;
	}
	else if (H > 13000.0 && H <= 26000.0)
	{
		DensutyAir = log(-298.95e-12*H*H + (12.780e+3) / H + 28.677e-6*H);
	}
	else if (H > 26000.0 && H <= 39000.0)
	{
		DensutyAir = exp(-7.9310*log(H) - 94590 / H + 80.891);
	}
	else if (H > 39000.0 && H <= 49000.0)
	{
		DensutyAir = exp(1.4632e-9*H*H + 264.69e-3*log(H) - 266.74e-6*H);
	}
	else
	{
		DensutyAir = 0.0;
	}
	return DensutyAir;
};
//������� ���������� ��� ��������� (������� ����������� �� ����� �����, �����, ������)
double cx(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	//����
	double m[10] = { 0.100, 0.500, 0.800, 0.900, 1.100, 1.200, 1.500, 1.900, 3.000, 100.000 };
	double al[5] = { 0.000, 2.000, 3.100, 5.000, 9.000 };
	double c[5][10] = { 0.521, 0.677, 0.747, 0.759, 0.824, 0.819, 0.803, 0.714, 0.549, 0.453,
						0.634, 0.678, 0.748, 0.759, 0.824, 0.820, 0.818, 0.727, 0.558, 0.460,
						0.637, 0.683, 0.753, 0.765, 0.826, 0.821, 0.850, 0.756, 0.578, 0.475,
						0.645, 0.693, 0.765, 0.776, 0.828, 0.824, 0.881, 0.783, 0.597, 0.491,
						0.672, 0.731, 0.806, 0.818, 0.837, 0.833, 0.950, 0.844, 0.642, 0.527 };

	alfa_ = alfa_ * 180 / pi;
	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	} while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
};

double cy_a(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { 1.250, 1.310, 1.418, 1.335, 1.266, 1.200, 1.150, 1.103, 0.979, 0.756,
						1.262, 1.331, 1.412, 1.329, 1.261, 1.196, 1.148, 1.102, 0.975, 0.746,
						1.289, 1.358, 1.412, 1.329, 1.261, 1.196, 1.148, 1.102, 0.975, 0.746,
						1.315, 1.386, 1.412, 1.329, 1.261, 1.196, 1.148, 1.102, 0.975, 0.746,
						1.447, 1.524, 1.412, 1.329, 1.261, 1.196, 1.148, 1.102, 0.975, 0.746 };

	alfa_ = alfa_ * 180 / pi;

	if (alfa_<0)
		alfa_ = -alfa_;

	if (alfa_ >= 8.0)
		alfa_ = 7.999;
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
};

double cz_b(double M, double alfa_)
{
	return -cy_a(M, alfa_);
};

double cy_d(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { 0.116, 0.117, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119,
						0.118, 0.119, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119,
						0.120, 0.121, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119,
						0.122, 0.123, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119,
						0.133, 0.134, 0.118, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119, 0.119 };

	alfa_ = alfa_ * 180 / pi;

	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef*1.1;
};

double cz_d(double M, double alfa_)
{
	return -cy_d(M, alfa_);
};
double mx_wx(double M, double alfa_)
{
	return -0.005*0.6786;
};
double mz_wz(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1, zn = 1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { -3.281, -3.035, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399,
						-3.357, -3.106, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399,
						-3.432, -3.177, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399,
						-3.508, -3.249, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399,
						-3.888, -3.605, -3.480, -3.789, -3.727, -3.525, -3.398, -3.281, -2.990, -2.399 };
	alfa_ = alfa_ * 180 / pi;

	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
};
double mz_a(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { -0.099, -0.084, -0.095, -0.104, -0.101, -0.096, -0.092, -0.089, -0.078, -0.054,
						-0.098, -0.084, -0.094, -0.103, -0.100, -0.095, -0.092, -0.089, -0.077, -0.052,
						-0.099, -0.085, -0.094, -0.103, -0.100, -0.095, -0.092, -0.089, -0.077, -0.052,
						-0.101, -0.086, -0.094, -0.103, -0.100, -0.095, -0.092, -0.089, -0.077, -0.052,
						-0.108, -0.091, -0.094, -0.103, -0.100, -0.095, -0.092, -0.089, -0.077, -0.052 };

	alfa_ = alfa_ * 180 / pi;

	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
};

double mz_d(double M, double alfa_)
{
	double v, n, koef;
	int j = -1, i = -1;
	double m[10] = { 0.500, 0.900, 1.200, 1.500, 1.800, 2.100, 2.300, 2.500, 3.000, 100.000 };
	double al[5] = { 0.000, 1.000, 2.000, 3.000, 8.000 };
	double c[5][10] = { 0.081, 0.083, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095,
						0.091, 0.094, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095,
						0.102, 0.104, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095,
						0.112, 0.115, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095,
						0.163, 0.168, 0.090, 0.093, 0.095, 0.097, 0.097, 0.097, 0.097, 0.095 };

	alfa_ = alfa_ * 180 / pi;

	if (alfa_ < 0)
	{
		alfa_ = -alfa_;
	}
	if (alfa_ >= 8.0)
	{
		alfa_ = 7.999;
	}
	do
	{
		j++;
	}
	while (!(M >= m[j] && M<m[j + 1]));
	do
	{
		i++;
	}
	while (!(alfa_ >= al[i] && alfa_<al[i + 1]));
	n = (c[i][j + 1] - c[i][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i][j];
	v = (c[i + 1][j + 1] - c[i + 1][j])*(M - m[j]) / (m[j + 1] - m[j]) + c[i + 1][j];
	koef = (v - n)*(alfa_ - al[i]) / (al[i + 1] - al[i]) + n;
	return koef;
}; 
double my_wy(double M, double alfa_)
{
	return mz_wz(M, alfa_);
};
double my_b(double M, double alfa_)
{
	return mz_a(M, alfa_);
};
double my_d(double M, double alfa_)
{
	return mz_d(M, alfa_);
};
double mx_d(double M, double alfa_)
{
	return 0.1*mz_d(M, alfa_);
};
//������� ���������� ������� ����� ����������� � ������������ ���������
double viz_fi(double Y, double r_viz, double Yc)
{
	return asin((Yc - Y) / r_viz);
};

//������� ���������� ������� ����� ����������� � �������������� ���������
double viz_hi(double X, double Xc, double Z, double Zc)
{
	return -atan2(Zc - Z, Xc - X);
};

//������� ���������� ������� ����� �����������
double r_vizV(double X, double Y, double Z, double Xc, double Yc, double Zc)
{
	return sqrt(pow((Yc - Y), 2) + pow((Xc - X), 2) + pow((Zc - Z), 2));
};