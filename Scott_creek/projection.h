#include <iostream>
#include "res.h"
class Projection
{
public:
		void LL_2_UTM(double latitude, double longitude, int zone, double &North, double &East);
		void UTM_2_LL(double North, double East, int zone, double &latitude, double &longitude);
};

void Projection::LL_2_UTM(double latitude, double longitude, int zone, double &North, double &East)
{
	double long_origin = -183 + zone*6;
	double a = 6378137.0;
	double f = 1.0/298.257222101;
	double e2 = 2*f - f*f;
	double e = sqrt(e2);
	double ep2 = e2/(1-e2); 

	double N = a/(sqrt(1-e2*sin(latitude)*sin(latitude)));
	double T = tan(latitude)*tan(latitude);
	double C = ep2*cos(latitude)*cos(latitude);
	double A = cos(latitude)*(longitude-(long_origin*PI/180));

	//Compute Meridian Arc Length
	double term1 = 1 - e2/4 - (3*e2*e2)/64 - (5*e2*e2*e2)/256;
	double term2 = (3*e2)/8 + (3*e2*e2)/32 + (45*e2*e2*e2)/1024;
	double term3 = (15*e2*e2)/256 + (45*e2*e2*e2)/1024;
	double term4 = (35*e2*e2*e2)/3072;

	double M = a*(term1*latitude - term2*sin(2*latitude) + term3*sin(4*latitude) - term4*sin(6*latitude));

	double x1 = ((1-T+C)*A*A*A)/6;
	double x2 = ((5 - 18*T + T*T + 72*C - 58*ep2)*A*A*A*A*A)/120;
	double x = 0.9996*N*(A + x1 + x2);

	double y1 = (5 - T + 9*C + 4*C*C)*(A*A*A*A)/24;
	double y2 = (61 - 58*T + T*T + 600*C - 330*ep2)*(A*A*A*A*A*A)/720;
	double y3 = (A*A)/2 + y1 + y2;
	double y = 0.9996*(M + N*tan(latitude)*y3);

	North = y;
	East = x + 500000.0;

	return;
};

void Projection::UTM_2_LL(double North, double East, int zone, double &latitude, double &longitude)
{
	double k_0=0.9996;
	double E_0=500000.0;
	double Ep=East-E_0;
	double N_0=0.0;
	double lamda_0=(-183.0+double(zone*6))*PI/180.0;
	double phi_0=0.0;
	double a=6378137.0;
	double f=1.0/298.2572221012;
	double e2=2.0*f-f*f;
	double ep2=e2/(1.0-e2);
	double n=f/(2.0-f);
	double n2=n*n;
	double n3=n*n*n;
	double n4=n*n*n*n;
	double r=a*(1.0-n)*(1.0-n2)*(1.0+9.0/4.0*n2+225.0/64.0*n4);
	double u2=-3.0/2.0*n+9.0/16.0*n3;
	double u4=15.0*n2/16.0-15.0*n4/32.0;
	double u6=-35.0*n3/48.0;
	double u8=315.0*n4/512.0;
	double U0=2.0*(u2-2.0*u4+3.0*u6-4.0*u8);
	double U2=8.0*(u4-4.0*u6+10.0*u8);
	double U4=32.0*(u6-6.0*u8);
	double U6=128.0*u8;
	double v2=3.0*n/2.0-27.0*n3/32.0;
	double v4=21.0*n2/16.0-55.0*n4/32.0;
	double v6=151.0*n3/96.0;
	double v8=1097.0*n4/512.0;
	double V0=2.0*(v2-2.0*v4+3.0*v6-4.0*v8);
	double V2=8.0*(v4-4.0*v6+10.0*v8);
	double V4=32.0*(v6-6.0*v8);
	double V6=128.0*v8;
	double omega_0=phi_0+sin(phi_0)*cos(phi_0)*(U0+U2*cos(phi_0)*cos(phi_0)+U4*pow(cos(phi_0),4)+U6*pow(cos(phi_0),6));
	double S_0=k_0*omega_0*r;
	double omega=(North-N_0+S_0)/(k_0*r);
	double phi_f=omega+(sin(omega)*cos(omega))*(V0+V2*pow(cos(omega),2)+V4*pow(cos(omega),4)+V6*pow(cos(omega),6));
	double t_f=tan(phi_f);
	double ita_f2=ep2*pow(cos(phi_f),2);
	double R_f=k_0*a/sqrt(1.0-e2*pow(sin(phi_f),2));
	double Q=Ep/R_f;
	double t_f2=pow(t_f,2.0);
	double t_f4=pow(t_f,4.0);
	double t_f6=pow(t_f,6.0);
	double B2=-1.0/2.0*t_f*(1.0+ita_f2);
	double B4=-1.0/12.0*(5.0+3.0*t_f2+ita_f2*(1-9*t_f2)-4*pow(ita_f2,2.0));
	double B6=1/360*(61.0+90.0*t_f2+45*t_f4+ita_f2*(46-252*t_f2-90*t_f4));
	latitude=phi_f+B2*pow(Q,2.0)*(1+pow(Q,2.0)*(B4+B6*pow(Q,2.0)));
	double B3=-1.0/6.0*(1.0+2.0*pow(t_f,2)+ita_f2);
	double B5=1.0/120.0*(5.0+28.0*pow(t_f,2)+24.0*pow(t_f,4)+ita_f2*(6.0+8.0*t_f2));
	double B7=-1.0/5040.0*(61.0+662.0*t_f2+1320*t_f4+720*t_f6);
	double L=Q*(1+pow(Q,2)*(B3+pow(Q,2)*(B5+B7*pow(Q,2))));
	longitude=lamda_0+L/cos(phi_f);
	
	return;
};