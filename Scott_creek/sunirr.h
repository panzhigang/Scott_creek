#ifndef __SUNIRR_H
#define __SUNIRR_H
#endif

#include "res.h"
#include <iostream>

class Sunirr
{
public:
	void sunirr(double elevation_angle, double azimuth_angle, double elevation, double* normal, double *irr);
};

void Sunirr::sunirr(double elevation_angle, double azimuth_angle, double elevation, double* normal, double *irr)
{
	double AM=1.0/(cos(deg2rad*(90.0-elevation_angle))+0.50572*(pow(6.07995+elevation_angle,-1.6364)));
	double I_d=1353.0*((1-0.14*elevation/1000.0)*pow(0.7,pow(AM,0.678))+0.14*elevation/1000.0);
	if(azimuth_angle<=270.0 && azimuth_angle>=0.0)
		azimuth_angle=90-azimuth_angle;
	else if(azimuth_angle>270 && azimuth_angle<=360)
		azimuth_angle=450-azimuth_angle;
	azimuth_angle*=deg2rad;
	elevation_angle*=deg2rad;
	double x=cos(elevation_angle)*cos(azimuth_angle);
	double y=cos(elevation_angle)*sin(azimuth_angle);
	double z=sin(elevation_angle);
	double s=x*normal[0]+y*normal[1]+z*normal[2];
	if(s<0)
		s=0;
	*irr=I_d*s;
};
