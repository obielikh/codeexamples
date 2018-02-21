/*
 ============================================================================
 Name        : e42.c
 Author      : Oleksii Bielikh
 Version     :
 Copyright   :
 Description : Brownian motion simulation
 ============================================================================
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define PI 3.141592653589



int main()
{
    int repeats = 100;
    double n;
    double m = 9.59 *PI* pow(10,-15);  //mass
    double kb = 1.38 * pow(10,-23);  //Boltzmann constant
    double temp = 297;
    double timestep = 0.01;
    double T = 48.5*pow(10,-3);  //relaxation time
    int steps = 1000;
    //variables related to coordinates have 'x' in their name, whereas
    //the ones with 'v' are related to velocity
    double sumx,sumv,sumx2,sumv2;
    double sqmeanv;
    double sqmeanx;
    double G1,G2;
    double vth;
    double c0;
    int i;
    int k;
    int u;

    double **a = (double **)malloc(steps * sizeof(double *));
    double **v = (double **)malloc(steps * sizeof(double *));
    double **x = (double **)malloc(steps * sizeof(double *));
	for (i=0; i<steps; i++){
		 a[i] = (double *)malloc(repeats * sizeof(double));
		 v[i] = (double *)malloc(repeats * sizeof(double));
		 x[i] = (double *)malloc(repeats * sizeof(double));}
    double *meanv = malloc((steps) * sizeof (double));
    double *meanx = malloc((steps) * sizeof (double));
    double *varv = malloc((steps) * sizeof (double));
    double *varx = malloc((steps) * sizeof (double));


    srand(time(NULL));

    n = 1/T;
    c0 = exp(-n*timestep);
    vth=1000*sqrt((kb*temp)/m);

    for (u = 0; u < repeats; u++){
        a[0][u] = -36 * (PI * PI) * 0.1;
        v[0][u] = 2.0;
        //vn[0][u] = 0;
        x[0][u] = 0.1;

        //verlet algorithm
        for (i = 1; i < steps; i++)
		{
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;
			G1 = (sqrt(log(U1) * (-2))) * (cos(2 * PI * U2)) ;
			G2 = (sqrt(log(U1) * (-2))) * (cos(2 * PI * U2)) ;
			v[i][u] = 0.5 * a[i-1][u]*timestep + sqrt(c0) * v[i-1][u] + vth * sqrt(1-c0) * G1;
			x[i][u] = x[i-1][u] + v[i][u] * timestep;
			a[i][u] = -36 * (PI * PI) * x[i][u];
			v[i][u] = 0.5 * a[i][u] * timestep + sqrt(c0) * v[i][u] + vth * sqrt(1-c0) * G2;
		}
	}
  //calculating mean values
	for(i =0; i < steps; i++)
	{
		sumv=0;
		sumv2=0;
		sumx=0;
		sumx2=0;
		for (u = 0; u < repeats; u++){
			sumv += v[i][u];
			sumv2 += (v[i][u]*v[i][u]);
			sumx += x[i][u];
			sumx2 += (x[i][u]*x[i][u]);
		}
		meanv[i] = sumv / repeats;
		sqmeanv = sumv2 / repeats;
		varv[i] = sqrt((sqmeanv - (meanv[i]*meanv[i])) / repeats);
		meanx[i] = sumx / repeats;
		sqmeanx = sumx2 / repeats;
		varx[i] = sqrt((sqmeanx - (meanx[i]*meanx[i])) / repeats);
	}

  //printing and cleaning

	FILE *file, *file2, *file3, *file4,*file5,*file6;
	file = fopen("ex11.dat", "w");
	file2 = fopen("ex22.dat", "w");
	file3=fopen("ex111.dat", "w");
	file4 = fopen("ex222.dat", "w");
	file5=fopen("ex1.dat", "w");
	file6 = fopen("ex2.dat", "w");
	for (k = 0; k < steps; k++) {
		fprintf(file, "%e\t%e\t%e\n",  (k*timestep), (meanv[k]), (varv[k]));
		fprintf(file2, "%e\t%e\t%e\n",  (k*timestep), (meanx[k]), (varx[k]));

		fprintf(file5, "%e\t%e\t%e\t%e\t%e\t%e\n",  (k*timestep), (v[k][1]),(v[k][2]),(v[k][3]),(v[k][4]),(v[k][5]));
		fprintf(file6, "%e\t%e\t%e\t%e\t%e\t%e\n",  (k*timestep), (x[k][1]),(x[k][2]),(x[k][3]),(x[k][4]),(x[k][5]));
	}
	fclose(file);
	fclose(file2);
	fclose(file5);
	fclose(file6);
	for (k=0; k<repeats; k++){
		fprintf(file3, "%e\t%e\n",  (x[485][k]), (v[485][k]));
		fprintf(file4, "%e\t%e\n",  (x[970][k]), (v[970][k]));
	}
	printf("done");
	fclose(file3);
	fclose(file4);

   	free(meanv); meanv = NULL;
	free(meanx); meanx = NULL;
	free(varv); varv = NULL;
	free(varx); varx = NULL;

	for(i = 0; i < steps; i++){
	free(a[i]);
	  a[i] = NULL;
	free(v[i]);
	  v[i] = NULL;
	free(x[i]);
	  x[i] = NULL;
	}
	free(a);
	  a= NULL;
	free(v);
	  v = NULL;
	free(x);
	  x = NULL;

	return 0;
}
