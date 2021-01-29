/*
 * LEAPINT.C: program to integrate hamiltonian system using leapfrog.
 */

#include <math.h>
#include <stdio.h>
#define MAXPNT  100				/* maximum number of points */

double G = 6.674*pow(10, -11);
double M = 1.9891*pow(10, 30);

double getVel();

void leapstep();				/* routine to take one step */

void accel();					/* accel. for harmonic osc. */

void printstate();				/* print out system state   */

void main(argc, argv)
int argc;
char *argv[];
{
    int n, mstep, nout, nstep;
    double x[MAXPNT][3], v[MAXPNT][3], tnow, dt;

    /* first, set up initial conditions */

    n = 9;

    int i, j;
    for(i = 0; i < n; i++){ //fill with zeros so I have to type less
        for (j = 0; j < 3; j++){
            x[i][j] = 0;
            v[i][j] = 0;
        }
    }
    				/* set number of points     */
    x[0][0] = 46.0 * pow(10, 9); // m
    x[1][0] = 107.5 * pow(10, 9);
    x[2][0] = 147.1 * pow(10, 9);
    x[3][0] = 206.6 * pow(10, 9);
    x[4][0] = 740.5 * pow(10, 9);
    x[5][0] = 1352.6 * pow(10, 9);
    x[6][0] = 2741.3 * pow(10, 9);
    x[7][0] = 4444.5 * pow(10, 9);
    x[8][0] = 4436.8 * pow(10, 9);

    v[0][1] = getVel( x[0][0], 46 * pow(10, 9), 69.8 * pow(10,9) );
    v[1][1] = getVel( x[1][0], 107.5 * pow(10, 9), 108.9 * pow(10,9) );
    v[2][1] = getVel( x[2][0], 147.1 * pow(10, 9), 152.1 * pow(10,9) );
    v[3][1] = getVel( x[3][0], 206.6 * pow(10, 9), 249.2 * pow(10,9) );
    v[4][1] = getVel( x[4][0], 740.5 * pow(10, 9), 816.6 * pow(10,9) );
    v[5][1] = getVel( x[5][0], 1352.6 * pow(10, 9), 1514.5 * pow(10,9) );
    v[6][1] = getVel( x[6][0], 2741.3 * pow(10, 9), 3003.6 * pow(10,9) );
    v[7][1] = getVel( x[7][0], 4444.5 * pow(10, 9), 4545.7 * pow(10,9) );
    v[8][1] = getVel( x[8][0], 4436.8 * pow(10, 9), 7375.9 * pow(10,9) );

    tnow = 0.0;					/* set initial time         */

    /* next, set integration parameters */

    mstep = 2560000;				/* number of steps to take  */
    nout = 20000;					/* steps between outputs    */
    dt = 3200;				/* timestep for integration */

    /* now, loop performing integration */

    for (nstep = 0; nstep < mstep; nstep++) {	/* loop mstep times in all  */
        if (nstep % nout == 0)			    /* if time to output state  */
            printstate(x, v, n, tnow);		/* then call output routine */
        leapstep(x, v, n, dt);			    /* take integration step    */
        tnow = tnow + dt;			        /* and update value of time */
    }   
    if (mstep % nout == 0)			    /* if last output wanted    */
	    printstate(x, v, n, tnow);		/* then output last step    */
}

/* Returns the instantaneous velocity for an object in orbit around the sun */
double getVel(d, peri, aphe)
double d;
double peri;
double aphe;
{
    return sqrt(G * M * ( 2/d - 2/(peri + aphe) ) );
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */
void leapstep(x, v, n, dt)
double x[MAXPNT][3];					/* positions of all points  */
double v[MAXPNT][3];					/* velocities of all points */
int n;						    /* number of points         */
double dt;					    /* timestep for integration */
{
    int i;
    double a[MAXPNT][3];

    accel(a, x, n);				        /* call acceleration code   */

    for (i = 0; i < n; i++){		        /* loop over all points...  */
        for(int j = 0; j < 3; j++)
	        v[i][j] = v[i][j] + 0.5 * dt * a[i][j];		/* advance vel by half-step */
    }
    for (i = 0; i < n; i++)	{		    /* loop over points again...*/
        for(int j = 0; j < 3; j++)
	        x[i][j] = x[i][j] + dt * v[i][j];		/* advance pos by full-step */
    }

    accel(a, x, n);				        /* call acceleration code   */
    
    for (i = 0; i < n; i++){			/* loop over all points...  */
        for (int j = 0; j < 3; j++)
	        v[i][j] = v[i][j] + 0.5 * dt * a[i][j];  /* and complete vel. step   */
    }		
}

/*
 * ACCEL: compute accelerations
 */

void accel(a, x, n)
double a[MAXPNT][3];					/* accelerations of points  */
double x[MAXPNT][3];					/* positions of points      */
int n;						/* number of points         */
{
    int i;
    int j;
    double mag;
    for (i = 0; i < n; i++){
        mag = sqrt( x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2] );
        for(j = 0; j < 3; j++){
	        a[i][j] = G * M * ( -x[i][j] / (mag*mag*mag) );				/* Making Newton proud? */
        }
    }
}

/*
 * PRINTSTATE: output system state variables.
 */

void printstate(x, v, n, tnow)
double x[MAXPNT][3];					/* positions of all points  */
double v[MAXPNT][3];				/* velocities of all points */
int n;						    /* number of points         */
double tnow;					/* current value of time    */
{
    int i;
    double r;
    double theta;
    double vr;
    double vtheta;

    for (i = 0; i < n; i++)			/* loop over all points...  */
    {
        r = sqrt( x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2] );
        theta = atan(x[i][1] / x[i][0]) + ((x[i][0]<0)?(M_PI):(0));
	    printf("%8.4f %4d %12.6f %12.6f\n",
            tnow, i, r, theta);
            // tnow, i, x[i][0], x[i][1], x[i][2], v[i][0], v[i][1], v[i][2]); /* time, index, x y z, vx vy vz */

    }
}
