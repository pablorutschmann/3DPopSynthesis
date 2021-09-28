#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <sys/stat.h> 
#include <stdlib.h>

using namespace std;

int fgfull(double *x, double *y, double *z, double *vx, double *vy, double *vz, double dt, double mu){

	/*-- KEPLER INTEGRATOR --*/

	/*
	INPUTS:
		- pointers to coordinates x, y, z, and velocities vx, vy, vz;
		- timestep
		- mu = G_grav * M_planet

	OUTPUTS:
		- return 0 if the procedure succeeds, 1 if it does not (satellite ejected)

	OTHER:
	    - the function assigns new values to coordinates and velocities
	*/

	double xi, yi, zi, vxi, vyi, vzi;
	double f,g,fd,gd;                               /* Gauss's f, g, fdot and gdot */
	double rsq,vsq,ir;
	double u;                                       /* r v cos(phi) */
	double ia;                                      //inverse of a
	double ria;
	double air_;
	double e;                                        /* eccentricity */
	double ec,es;                             	  /* e cos(E), e sin(E) */
	double ien;                                   /* inverse mean motion */
	double en;                                    /* mean motion */
	double dec;                                      /* delta E */
	double dm;
	double mw;                                       /* minus function to zero */
	double wp;                                       /* first derivative */
	double iwp;
	double wpp;                                      /* second derivative */
	double wppp;                                     /* third derivative */
	double dx,s,c;
	double t1;
	double tx, ty, tz;
	const double DOUBLE_EPS = 1.2e-16;
	double converge;
	double UP = 2*M_PI;
	double LOW = -2*M_PI;
	double next;
	int i;
	/*
	* Evaluate some orbital quantites.
	*/

	xi = *x;
	yi = *y;
	zi = *z;
	vxi = *vx;
	vyi = *vy;
	vzi = *vz;

	rsq = xi*xi + yi*yi + zi*zi;
	vsq = vxi*vxi + vyi*vyi + vzi*vzi;
	u =  xi*vxi + yi*vyi + zi*vzi;
	ir = 1.0 / sqrt(rsq);
	ia = 2.0*ir-vsq/mu;
	
	if(ia > 0)
	{
		t1 = ia*ia;
		ria = rsq*ir*ia;
		ien = 1.0 / sqrt(mu*t1*ia);
		en = 1.0/ien;
		ec = 1.0-ria;
		es = u*t1*ien;
		e = sqrt(ec*ec + es*es);
		dm = en * dt - es;
		if ((es*cos(dm)+ec*sin(dm)) > 0){
			dec = dm + 0.85*e;
		}
		else dec = dm - 0.85*e;
		converge = fabs(en * dt *DOUBLE_EPS);
		for(i = 0; i < 128; ++i) {

			s = sin(dec);
			c = cos(dec);
			//sincos(dec, &s, &c);
			wpp = ec*s + es*c;
			wppp = ec*c - es*s;
			mw = dm - dec + wpp;
			if(mw < 0){
				UP = dec;
			}
			else LOW = dec;
			wp = 1.0 - wppp;
			wpp *= 0.5;
			dx = mw/wp;
			dx = mw/(wp + dx*wpp);
			dx = mw/(wp + dx*(wpp + (1.0/6.0)*dx*wppp));
			next = dec + dx;
			if (fabs(dx) <= converge) break;
			if(next > LOW && next < UP){
				dec = next;
			}
			else dec = 0.5*(LOW + UP);
			if (dec==LOW || dec==UP) break;
		}
		/*
		if(i < 127){
			iwp = 1.0/wp;
			air_ = -1.0/ria;
			t1 = (1.0-c);
			f = 1.0 + air_*t1;
			g = dt + (s-dec)*ien;
			fd = air_*iwp*s*en;
			gd = 1.0 - iwp*t1;

			tx = f*xi+g*vxi;
			ty = f*yi+g*vyi;
			tz = f*zi+g*vzi;
			vxi = fd*xi+gd*vxi;
			vyi = fd*yi+gd*vyi;
			vzi = fd*zi+gd*vzi;
			xi = tx;
			yi = ty;
			zi = tz;

			*x = xi;
			*y = yi;
			*z = zi;
			*vx = vxi;
			*vy = vyi;
			*vz = vzi;
		}
		else{
            cout << "not converged\n";
		}
		*/

		iwp = 1.0/wp;
		air_ = -1.0/ria;
		t1 = (1.0-c);
		f = 1.0 + air_*t1;
		g = dt + (s-dec)*ien;
		fd = air_*iwp*s*en;
		gd = 1.0 - iwp*t1;

		tx = f*xi+g*vxi;
		ty = f*yi+g*vyi;
		tz = f*zi+g*vzi;
		vxi = fd*xi+gd*vxi;
		vyi = fd*yi+gd*vyi;
		vzi = fd*zi+gd*vzi;
		xi = tx;
		yi = ty;
		zi = tz;

		*x = xi;
		*y = yi;
		*z = zi;
		*vx = vxi;
		*vy = vyi;
		*vz = vzi;
		
		return 0;
	}
	else
	{
        return 1;
	}
}