#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "csgraphics.h"

#define MAXR 50
#define MAXACTIVE 10000
#define TWOPI 6.283185307179586476925286766559
#define S 0.333333333333333333333333333333
#define P 0.18
#define V 1
#define DT 10

void parameter();
void initial();
void evolve();
void create(int r, int a);
void plot_spiral();

double t;
int active_r[MAXACTIVE];
int active_a[MAXACTIVE];
int cell[MAXR][MAXR*6+1];
int nring, nactive;
double s, p, v, dt, two_pi;

int main()
{
    srand(time(NULL));
    parameter();
    initial();
    while (!kbhit()) {
        evolve();
        plot_spiral();
        t += dt;
    }
    return 0;
}

void parameter()
{
    nring = 50;
    nactive = 200;
    v = V;
    p = P;
    dt = DT;
    s = S;
    two_pi = TWOPI;
}

void initial()
{
	int i, x, y, r, a;
	int active_a[10000], active_r[10000];
	int nactive = 10000;
	int cell[nring][6*nring];

	randomize();
	// randomly activate cells
	i = 0;
	while (i < nactive) {
		do {
			x = int(nring*rnd) + 1;
			y = int(nring*rnd) + 1;
			r = int(sqr(x*x + y*y)) + 1;
		} while (r <= nring);
		a = int(6*r*rnd) + 1;		// array index corresponds to an angle
		if (cell[r][a] == 0) {
			i = i + 1;
			active_a[i] = a;		// location of active region
			active_r[i] = r;
			// activate region, stars live for 15 time steps 
			cell[r][a] = 15;
		}
	}
	t = 0;		// initial time
	compute_aspect_ratio(nring, xwin,ywin); 
	set_window(-xwin, xwin, -ywin, ywin);
	set_background_color("white");
	set_color("black");
}PROGRAM galaxy
! model proposed by Schulman and Seiden
LIBRARY "csgraphics" 
PUBLIC t,active_r(10000), active_a(10000), cell(0 to 50,300)
DECLARE PUBLIC dt
CALL parameter
CALL initial

DO
	CALL evolve
	CALL plot spiral 
	LET t = t + dt
LOOP UNTIL key input
END

SUB parameter
	PUBLIC nring.nactive,v,p,dt,s,two_pi
	LET nring = 50			! number of rings 
	LET nactive = 200		! number of initial active cells
	LET v = 1			! circular velocity 
	LET p = 0.18			! star formation probability		
	LET dt = 10			! time step
	LET s = 2*pi/6			! cell width
	LET two_pi = 2*pi

END SUB

SUB initial
	DECLARE PUBLIC t, active_r(), active_a(), cell(,), nring, nactive
	RANDOMIZE
	MAT cell = 0
	! randomly activate cells
	LET i = 0
	DO
		DO
			LET x = int(nring*rnd) + 1
			LET y = int(nring*rnd) + 1
			LET r = int(sqr(x*x + y*y)) + 1
		LOOP until r <= nring
		LET a = int(6*r*rnd) + 1		! array index corresponds to an angle
		IF cell (r,a) = 0 then 
			LET i = i + 1
			LET active_a(i) = a		! location of active region
			LET active_r(i) = r
			! activate region, stars live for 15 time steps 
			LET cell (r,a) = 15
		END IF
	LOOP until i = nactive
	LET t = 0		! initial time
	CALL compute_aspect_ratio(nring, xwin,ywin) 
	SET WINDOW -xwin, xwin, -ywin, ywin
	SET BACKGROUND COLOR "white"
	SET COLOR "black"
END SUB

SUB evolve
	DECLARE PUBLIC t, active_r(), active_a(), cell (,) 
	DECLARE PUBLIC nring, nactive,v,p,s,two_pi
	DIM newactive_r(10000), newactive_a(10000)
	! number of active star clusters for next time step
	LET newactive = 0 
	FOR i = 1 to nactive
		LET r = active_r(i)
		LET a = active_a(i)
		! activate neighboring cells of same ring
		CALL create (r, a+1, newactive, newactive_r(), newactive_a()) 
		CALL create (r, a-1, newactive, newactive_r(), newactive_a())
		LET angle mod ((a*s + v*t)/r, two_pi)
		IF r < nring then		! activate cells in next larger ring
			LET wt = mod(v*t/(r+1), two_pi)
			LET ap = int((angle - wt)*(r+1)/s) 
			LET ap = mod(ap, 6*(r+1))
			CALL create (r+1, ap, newactive, newactive_r(), newactive_a())
			CALL create (r+1, ap+1,newactive, newactive_r(), newactive_a())
		END IF
		IF r > 1 then		! activate cells in next smaller ring
			LET wt = mod (v*t/(r-1), two_pi)
			LET am = int((angle - wt)*(r-1)/s) 
			LET am = mod (am, 6*(r-1))
			CALL create (r-1, am, newactive, newactive_r(),newactive_a()) 
			CALL create (r-1, an+1,newactive, newactive_r(), newactive_a())
		END IF
	NEXT i
	LET nactive = newactive 
	MAT active_r= newactive_r
	MAT active a newactive_a
END SUB

SUB create (r,a,nevactive, newactive_r(),newactive_a()) 
	! create star clusters
	DECLARE PUBLIC p. cell(,)
	IF a < 1 then LET a = a + 6*r
	IF rnd < p and cell (r,a) <> 15 then
		LET newactive = newactive + 1
		LET newactive_a(newactive) = a
		LET newactive_r(newactive) = r 
		LET cell (r,a) = 15	! activate cell
	END IF
END SUB

SUB plot spiral
	DECLARE PUBLIC t,nring, nactive,v,s, cell(,) 
	CLEAR
	SET CURSOR 1,1
	PRINT "number of active star clusters ="; nactive
	FOR r = 1 to nring
		FOR a = 1 to 6*r
			IF cell (r,a) > 0 then
				LET theta = (a*s + v*t)/r 
				LET x = r*cos(theta)
				LET y = r*sin(theta)
				LET plotsize cell (r,a)/30 
				BOX AREA x-plotsize,x+plotsize,y-plotsize,y+plotsize
				! reduce star clusters lifetime 
				LET cell(r,a) = cell (r,a) - 1
			END IF
		NEXT a
	NEXT r
END SUB