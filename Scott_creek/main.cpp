#include "projection.h"
#include "sunpos.h"
#include "sunirr.h"
void main()
{
	//Projection pro;
	//double Lat, Lon;
	//pro.UTM_2_LL(3320469.28642977,307084.894817944,15,Lat,Lon);
	cTime time;
	cLocation location;
	cSunCoordinates sun;
	time.iYear=2014;
	time.iMonth=3;
	time.iDay=20;
	time.dHours=20.0;
	time.dMinutes=0.0;
	time.dSeconds=0.0;
	location.dLatitude=30.0;
	location.dLongitude=-95.0;
	Sunpos sunp;
	sunp.sunpos(time,location,&sun);
	Sunirr irr;
	//double ele=60.0;
	//double azi=90.0;
	double normal[3];
	normal[0]=0.0;
	normal[1]=0.0;
	normal[2]=1.0;
	double irradiation;
	irr.sunirr(90.0-sun.dZenithAngle,sun.dAzimuth,0.0,normal,&irradiation);
}