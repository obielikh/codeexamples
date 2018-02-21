/*
 ============================================================================
 Name        : E1T1.c
 Author      : Oleksii Bielikh
 Version     :
 Copyright   :
 Description : Harmonic oscillators
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

//Signal generation
void func1( int n, double f, double a, double fi, double delta, char name[])
{
	double val = 0;
	double t[n]; //time
	double h[n];  //signal
	for( int i = 0; i < n; i++)
	{
		t[i] = val;
		val += delta;
	}

	for(int j = 0; j < n; j++)
	{
		h[j] = a*cos(2*M_PI*f*t[j] + fi);
	}
	FILE *fil = fopen(name, "w");
	for (int k = 0; k < n; k++) {
	  fprintf(fil, "%.3g\t%.5g\n", t[k], h[k]);
	}
	fclose(fil);
}
//Generation of corresponding power spectrum
void func2( int n, double f, double a, double fi, double delta, char name[])
{
	double val = 0;
	double t[n];
	double h[n];
	float y[n];  //absolute value
	float p[n];  //power
	double deltaf;
	float hre[n], him[n];  //real and imaginary part
	deltaf = 1/(n*delta);
	float df[n];  //frequency
	for( int i = 0; i < n; i++)
	{
		t[i] = val;
		val += delta;
	}

	for(int j = 0; j < n; j++)
	{
		h[j] = a*cos(2*M_PI*f*t[j] + fi);
	}

	 for (int ky = 0 ; ky <n ; ++ky)
	{
		hre[ky] = 0;
		for (int q = 0 ; q < n ; ++q)
		{
			hre[ky] += h[q] * cos(q * ky * (2*M_PI) / n);
		}
		him[ky] = 0;
		for (int q = 0 ; q < n ; ++q)
		{
			him[ky] -= h[q] * sin(q * ky * (2*M_PI) / n);
		}
		y[ky] = hre[ky]*hre[ky] + him[ky]*him[ky];
	}
	 for (int i = 0; i < n; i++)
	{
		p[i] = y[i]/n;
	}
	 val = 0;
	for (int i = 0; i < n; i++)
	{
	 	df[i] = val;
		val += deltaf;
	}
	FILE *fil = fopen(name, "w");
	for (int k = 0; k < n; k++) {
	  fprintf(fil, "%g\t%f\n", df[k], p[k]);
	}
	fclose(fil);
}
//Making power spectrum plot more convenient
void func3( int n, double f, double a, double fi, double delta, char name[],char name2[])
{
	double val = 0;
	double t[n];
	double h[n];
	float y[n];
	float p[n];
	double deltaf;
	float hre[n], him[n];
	deltaf = 1/(n*delta);
	float df[n];
	for( int i = 0; i < n; i++)
	{
		t[i] = val;
		val += delta;
	}

	for(int j = 0; j < n; j++)
	{
		h[j] = a*cos(2*M_PI*f*t[j] + fi);
	}

	 for (int ky = 0 ; ky <n ; ++ky)
	{
		hre[ky] = 0;
		for (int q = 0 ; q < n ; ++q)
		{
			hre[ky] += h[q] * cos(q * ky * (2*M_PI) / n);
		}
		him[ky] = 0;
		for (int q = 0 ; q < n ; ++q)
		{
			him[ky] -= h[q] * sin(q * ky * (2*M_PI) / n);
		}
		y[ky] = hre[ky]*hre[ky] + him[ky]*him[ky];
	}
	 for (int i = 0; i < n; i++)
	{
		p[i] = y[i]/n;
	}
	 val = 0;
	for (int i = 0; i < n; i++)
	{
	 	df[i] = val;
		val += deltaf;
	}
	double t2[n];
	double h2[n];
	float p2[n];
	float df2[n];
	int in = 0;
	for(int i = n/2; i<n; i++)
	{
		t2[in] = t[i];
		h2[in] = h[i];
		p2[in] = p[i];
		df2[in] = df[i];
		in++;
	}
	for(int i = 0; i < n/2; i++)
	{
		t2[in] = t[i];
		h2[in] = h[i];
		p2[in] = p[i];
		df2[in] = df[i];
		in++;
	}
	FILE *fil = fopen(name, "w");
	for (int k = 0; k < n; k++) {
	  fprintf(fil, "%g\t%f\n", t2[k], h2[k]);
	}
	fclose(fil);
	FILE *fil2 = fopen(name2, "w");
	for (int k = 0; k < n; k++) {
	  fprintf(fil2, "%g\t%f\n", df2[k], p2[k]);
	}
	fclose(fil2);
}
//Modification for new signal
void func4( int n, double f1, double a1, double fi1, double f2, double a2, double fi2,double delta, char name[])
{
	double val = 0;
	double t[n];
	double h[n];
	float y[n];
	float p[n];
	double deltaf;
	float hre[n], him[n];
	deltaf = 1/(n*delta);
	float df[n];
	for( int i = 0; i < n; i++)
	{
		t[i] = val;
		val += delta;
	}

	for(int j = 0; j < n; j++)
	{
		h[j] = a1*cos(2*M_PI*f1*t[j] + fi1) + a2*cos(2*M_PI*f2*t[j] + fi2);
	}

	 for (int ky = 0 ; ky <n ; ++ky)
	{
		hre[ky] = 0;
		for (int q = 0 ; q < n ; ++q)
		{
			hre[ky] += h[q] * cos(q * ky * (2*M_PI) / n);
		}
		him[ky] = 0;
		for (int q = 0 ; q < n ; ++q)
		{
			him[ky] -= h[q] * sin(q * ky * (2*M_PI) / n);
		}
		y[ky] = hre[ky]*hre[ky] + him[ky]*him[ky];
	}
	 for (int i = 0; i < n; i++)
	{
		p[i] = y[i]/n;
	}
	 val = 0;
	for (int i = 0; i < n; i++)
	{
		df[i] = val;
		val += deltaf;
	}
	double t2[n];
	double h2[n];
	float p2[n];
	float df2[n];
	int in = 0;
	for(int i = n/2; i<n; i++)
	{
		t2[in] = t[i];
		h2[in] = h[i];
		p2[in] = p[i];
		df2[in] = df[i];
		in++;
	}
	for(int i = 0; i < n/2; i++)
	{
		t2[in] = t[i];
		h2[in] = h[i];
		p2[in] = p[i];
		df2[in] = df[i];
		in++;
	}
	FILE *fil = fopen(name, "w");
	for (int k = 0; k < n; k++) {
	  fprintf(fil, "%g\t%f\n", df2[k], p2[k]);
	}
	fclose(fil);
}
int main(void) {
	func1(250, 2.0, 1.0, 0.0, 0.1, "1_1");
	func1(250, 1.0, 1.0, 0.0, 0.1, "1_2");
	func1(250, 1.0, 1.0, M_PI/2, 0.1, "1_3");

	func2(250, 2.0, 1.0, 0.0, 0.1, "2_1");
	func2(250, 1.0, 1.0, 0.0, 0.1, "2_2");
	func2(250, 1.0, 1.0, M_PI/2, 0.1, "2_3");

	func3(250, 2.0, 1.0, 0.0, 0.1, "3_1_1","3_1_2");
	func3(250, 1.0, 1.0, 0.0, 0.1, "3_2_1","3_2_2");
	func3(250, 1.0, 1.0, M_PI/2, 0.1, "3_3_1","3_3_2");

	func3(258, 2.0, 1.0, 0.0, 0.1, "3_4_1","3_4_2");
	func3(258, 1.0, 1.0, 0.0, 0.1, "3_5_1","3_5_2");
	func3(258, 1.0, 1.0, M_PI/2, 0.1, "3_6_1","3_6_2");

	func3(260, 2.0, 1.0, 0.0, 0.1, "3_7_1","3_7_2");
	func3(260, 1.0, 1.0, 0.0, 0.1, "3_8_1","3_8_2");
	func3(260, 1.0, 1.0, M_PI/2, 0.1, "3_9_1","3_9_2");

	func4(250, 2.0, 1.0, 0.0, 6.0, 1.0, 0.0, 0.1, "4_1");
	func4(250, 2.0, 1.0, 0.0, 6.0, 1.0, 0.0, 0.05, "4_2");
	return EXIT_SUCCESS;
}
