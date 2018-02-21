/*
 ============================================================================
 Name        : t1.c
 Author      : Oleksii Bielikh
 Version     :
 Copyright   :
 Description : Monte Carlo integration
 ============================================================================
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#define PI 3.141592653589

//one-dimensional integral
void t1()
{
	int case1 = 10;
	int case2 = 100;
	int case3 = 1000;
	int case4 = 10000;
	double ran, sum, sum2;
	double mean1, mean2, mean3, mean4;
	double sqmean1, sqmean2, sqmean3, sqmean4;
	double var1, var2, var3, var4;
	double case1Arr[case1],case2Arr[case2],case3Arr[case3],case4Arr[case4];
	double f1[case1], f2[case2], f3[case3], f4[case4];
	double u;
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));

	//random initialization
	for(int i = 0; i < case1; i++)
	{
	  ran = gsl_rng_uniform(q);
	  case1Arr[i] = ran;
	}
	for(int i = 0; i < case2; i++)
	{
		ran = gsl_rng_uniform(q);
	  case2Arr[i] = ran;
	}
	for(int i = 0; i < case3; i++)
	{
		ran = gsl_rng_uniform(q);
	  case3Arr[i] = ran;
	}
	for(int i = 0; i < case4; i++)
	{
		ran = gsl_rng_uniform(q);
	  case4Arr[i] = ran;
	}
	gsl_rng_free(q);
	//calculating the integral and mean values
	sum = 0;
	sum2 = 0;
	for(int i =0; i < case1; i++)
	{
	  f1[i] = case1Arr[i] * (1 - case1Arr[i]);
	  sum += f1[i];
	  sum2 += (f1[i] * f1[i]);
	}
	mean1 = sum / case1;
	sqmean1 = sum2 / case1;
	var1 = sqrt((sqmean1 - (mean1*mean1)) / case1);

	sum = 0;
	sum2 = 0;
	for(int i =0; i < case2; i++)
	{
	  f2[i] = case2Arr[i] * (1 - case2Arr[i]);
	  sum += f2[i];
	  sum2 += (f2[i] * f2[i]);
	}
	mean2 = sum / case2;
	sqmean2 = sum2 / case2;
	var2 = sqrt((sqmean2 - (mean2*mean2)) / case2);

	sum = 0;
	sum2 = 0;
	for(int i =0; i < case3; i++)
	{
	  f3[i] = case3Arr[i] * (1 - case3Arr[i]);
	  sum += f3[i];
	  sum2 += (f3[i] * f3[i]);
	}
	mean3 = sum / case3;
	sqmean3 = sum2 / case3;
	var3 = sqrt((sqmean3 - (mean3*mean3)) / case3);

	sum = 0;
	sum2 = 0;
	for(int i =0; i < case4; i++)
	{
	  f4[i] = case4Arr[i] * (1 - case4Arr[i]);
	  sum += f4[i];
	  sum2 += (f4[i] * f4[i]);
	}
	mean4 = sum / case4;
	sqmean4 = sum2 / case4;
	var4 = sqrt((sqmean4 - (mean4*mean4)) / case4);
	//printing
	FILE *file;
	file = fopen("e3t1.dat", "w");
	fprintf(file, "%e\t%e\n", mean1, var1);
	fprintf(file, "%e\t%e\n", mean2, var2);
	fprintf(file, "%e\t%e\n", mean3, var3);
	fprintf(file, "%e\t%e\n", mean4, var4);
	fclose(file);
	printf("done");
}
//importance sampling
void t2()
{
	int case1 = 10;
	int case2 = 100;
	int case3 = 1000;
	int case4 = 10000;
	double ran, sum, sum2;
	double mean1, mean2, mean3, mean4;
	double sqmean1, sqmean2, sqmean3, sqmean4;
	double var1, var2, var3, var4;
	double case1Arr[case1],case2Arr[case2],case3Arr[case3],case4Arr[case4];
	double f1[case1], f2[case2], f3[case3], f4[case4];
	double u;
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));

//random initialization with weight function - transformation method
	for(int i = 0; i < case1; i++)
	{
		ran = gsl_rng_uniform(q);
		case1Arr[i] = acos(1-2*ran) / PI;
		//case1Arr[i]=asin(2*ran/PI)/PI;
	}
	for(int i = 0; i < case2; i++)
	{
		ran = gsl_rng_uniform(q);
		case2Arr[i] = acos(1-2*ran) / PI;
		//case2Arr[i]=asin(2*ran/PI)/PI;
	}
	for(int i = 0; i < case3; i++)
	{
		ran = gsl_rng_uniform(q);
		case3Arr[i] = acos(1-2*ran) / PI;
		//case3Arr[i]=asin(2*ran/PI)/PI;
	}
	for(int i = 0; i < case4; i++)
	{
		ran = gsl_rng_uniform(q);
		case4Arr[i] = acos(1-2*ran) / PI;
		//case4Arr[i]=asin(2*ran/PI)/PI;
	}

	sum = 0;
	sum2 = 0;
	//calculating the integral and mean values
	for(int i =0; i < case1; i++)
	{
	  f1[i] = case1Arr[i] * (1 - case1Arr[i]);
	  f1[i]=f1[i]/((PI/2)*sin(PI*case1Arr[i]));
	  sum += f1[i];
	  sum2 += (f1[i] * f1[i]);
	}
	mean1 = sum / case1;
	sqmean1 = sum2 / case1;
	var1 = sqrt((sqmean1 - (mean1*mean1)) / case1);

	sum = 0;
	sum2 = 0;
	for(int i =0; i < case2; i++)
	{
	  f2[i] = case2Arr[i] * (1 - case2Arr[i]);
	  f2[i]=f2[i]/((PI/2)*sin(PI*case2Arr[i]));
	  sum += f2[i];
	  sum2 += (f2[i] * f2[i]);
	}
	mean2 = sum / case2;
	sqmean2 = sum2 / case2;
	var2 = sqrt((sqmean2 - (mean2*mean2)) / case2);

	sum = 0;
	sum2 = 0;
	for(int i =0; i < case3; i++)
	{
	  f3[i] = case3Arr[i] * (1 - case3Arr[i]);
	  f3[i]=f3[i]/((PI/2)*sin(PI*case3Arr[i]));
	  sum += f3[i];
	  sum2 += (f3[i] * f3[i]);
	}
	mean3 = sum / case3;
	sqmean3 = sum2 / case3;
	var3 = sqrt((sqmean3 - (mean3*mean3)) / case3);

	sum = 0;
	sum2 = 0;
	for(int i =0; i < case4; i++)
	{
	  f4[i] = case4Arr[i] * (1 - case4Arr[i]);
	  f4[i]=f4[i]/((PI/2)*sin(PI*case4Arr[i]));
	  sum += f4[i];
	  sum2 += (f4[i] * f4[i]);
	}
	mean4 = sum / case4;
	sqmean4 = sum2 / case4;
	var4 = sqrt((sqmean4 - (mean4*mean4)) / case4);

	FILE *file;
	file = fopen("e3t2.dat", "w");
	fprintf(file, "%e\t%e\n", mean1, var1);
	fprintf(file, "%e\t%e\n", mean2, var2);
	fprintf(file, "%e\t%e\n", mean3, var3);
	fprintf(file, "%e\t%e\n", mean4, var4);
	fclose(file);
	//printing
	int j;
	FILE *file2;
	file2 = fopen("e3t2_a1.dat", "w");
	for (j = 0; j < case1; j++)
	{
		fprintf(file2, "%e\n", case1Arr[j]);
	}
	fclose(file2);

	FILE *file3;
	file3 = fopen("e3t2_a2.dat", "w");
	for (j = 0; j < case2; j++)
	{
		fprintf(file3, "%e\n", case2Arr[j]);
	}
	fclose(file3);

	FILE *file4;
	file4 = fopen("e3t2_a3.dat", "w");
	for (j = 0; j < case3; j++)
	{
		fprintf(file4, "%e\n", case3Arr[j]);
	}
	fclose(file4);

	FILE *file5;
	file5 = fopen("e3t2_a4.dat", "w");
	for (j = 0; j < case4; j++)
	{
		fprintf(file5, "%e\n", case4Arr[j]);
	}
	fclose(file5);
	printf("done");
}
//three-dimensional integral
void t3()
{

	int count = 10000;
	double ran1, ran2, sum, mean, var;
	double x[count], y[count], z[count];
	double tempx, tempy, tempz;
	const gsl_rng_type *T;
	gsl_rng *q;
	double delta = 2.3;
	int no;
	int steps = 10000000;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));

	for(int i = 0; i < count; i++)
	{
		ran1 = gsl_rng_uniform(q);
		ran2 = gsl_rng_uniform(q);
		x[i] = ((sqrt(log(ran1) * (-2))) * (cos(2 * PI * ran2))) / sqrt(2);
		ran1 = gsl_rng_uniform(q);
		ran2 = gsl_rng_uniform(q);
		y[i] = ((sqrt(log(ran1) * (-2))) * (cos(2 * PI * ran2))) / sqrt(2);
		ran1 = gsl_rng_uniform(q);
		ran2 = gsl_rng_uniform(q);
		z[i] = ((sqrt(log(ran1) * (-2))) * (cos(2 * PI * ran2))) / sqrt(2);
	}


	FILE *file;
	file = fopen("e3t3.dat", "w");
	for (int i = 0; i < count; i++)
	{
		fprintf(file, "%e\t%e\t%e\n", x[i], y[i], z[i]);
	}
	fclose(file);
	//Metropolis
	no = 0;
	int index;
	for(int i = 0; i < steps; i++)
	{
		//Maybe use the same random value for x y and z updates
		ran1 = gsl_rng_uniform(q) * 1000;
		index = ran1;
		ran2 = gsl_rng_uniform(q);
		tempx = x[index] + (delta * (ran2 - 0.5));
		ran2 = gsl_rng_uniform(q);
		tempy = y[index] + (delta * (ran2 - 0.5));
		ran2 = gsl_rng_uniform(q);
		tempz = z[index] + (delta * (ran2 - 0.5));


		var = exp(-((tempx * tempx) + (tempy * tempy) + (tempz * tempz)) + ((x[index] * x[index]) + (y[index] * y[index]) + (z[index] * z[index])));
		ran2 = gsl_rng_uniform(q);
		if( var >= ran2)
		{
			no++;
			x[index] = tempx;
			y[index] = tempy;
			z[index] = tempz;
		}
	}
	printf("change %d" , no);

	FILE *file2;
	file = fopen("e3t3_2.dat", "w");
	for (int i = 0; i < count; i++)
	{
		fprintf(file, "%e\t%e\t%e\n", x[i], y[i], z[i]);
	}
	fclose(file);
	sum = 0;
	for(int i = 0; i < count; i++)
	{
		/*sum += (((x[i] * x[i]) + ((y[i] * y[i]) * (x[i] * x[i])) + ((z[i] * z[i]) * (y[i] * y[i]) * (x[i] * x[i]))) *
				exp(-((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i])))/
				((pow(PI , -1.5))*exp(-((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i])))));*/
		sum+=((x[i] * x[i]) + ((y[i] * y[i]) * (x[i] * x[i])) + ((z[i] * z[i]) * (y[i] * y[i]) * (x[i] * x[i])));


	}
	//mean = (sum / count) * (pow(PI , -1.5));
	mean = (sum / count);
	printf("\nMean: %e", mean);
}
//determining the statistical inefficiency
void t4()
{

	  int i, nbr_of_lines;
	  double mean1, mean2, sum, sum2;
	  FILE *in_file;
	  nbr_of_lines = 1e6; /* The number of lines in MC.txt. */
	  double *fik = malloc(((nbr_of_lines)) * sizeof (double));
	  double *data = malloc(((nbr_of_lines)) * sizeof (double));
	  sum = 0;
	  /* Read data from file. */
	  in_file = fopen("MC.txt","r");
	  for (i=0; i<nbr_of_lines; i++) {
			fscanf(in_file,"%lf",&data[i]);
			sum += data[i];
			sum2 += (data[i] * data[i]);
	  }
	  fclose(in_file);
	  mean1 = sum / nbr_of_lines;
	  mean2 = sum2 / nbr_of_lines;
	  printf("\nMean: %e", mean1);
	  printf("\nMean: %e", mean2);
	  /*FILE *ts;
		ts = fopen("test.dat", "w");*/
		//auto-correlation function
	  for(int k = 0; k < 20; k++)
	  {
		  sum = 0;
		  for (i=0; i<=nbr_of_lines; i++) {
			  if((i+k) > nbr_of_lines-1)
			  {
				  sum += (data[i] * data[ i + k - nbr_of_lines]);
			  }
			  else
			  {
				  sum += (data[i] * data[(i + k)]);
			  }
			  //fprintf(ts, "%e\n", sum);
		  }

		 fik[k] = ((sum / nbr_of_lines) - (mean1 * mean1)) / (mean2 - (mean1 * mean1));
	  }
	 // fclose(ts);

	FILE *file;
	file = fopen("e3t4_1.dat", "w");
	for (int j = 0; j < 20; j++)
	{
		fprintf(file, "%d\t%e\n", j, fik[j]);
	}
	fclose(file);

	//block averaging
	int b[16] = { 10, 16, 20, 25, 32, 40, 50, 64, 80, 100, 125, 160, 200, 250, 320, 400};
	double var1, var2, mean11, mean21;
	int temp;
	double s[16];
	for( int y = 0; y < 16; y++)
	{
		temp = b[y];
		double *arr = malloc(((nbr_of_lines / temp)) * sizeof (double));
		for(int j = 0; j < (nbr_of_lines / temp); j++)
		{
			sum =0;
			for( i = 0; i < temp; i++)
			{
				sum += data[i+(j*temp)];
			}
			arr[j] = sum / temp;
		}
		var1 = mean2 - (mean1*mean1);
		sum = 0;
		sum2 = 0;
		for(int p = 0; p < (nbr_of_lines / temp); p++)
		{
			sum += arr[p];
			sum2 += (arr[p] * arr[p]);
		}
		mean11 = sum / (nbr_of_lines /temp);
		mean21 = sum2 / (nbr_of_lines /temp);
		var2 = mean21 - (mean11*mean11);
		s[y] = (temp * var2)/var1;
		free(arr); arr = NULL;
	}

	printf("done");
	FILE *fileu;
	fileu = fopen("e3t4_2.dat", "w");
	for (int j = 0; j < 16; j++)
	{
		fprintf(fileu, "%d\t%e\n", b[j], s[j]);
	}
	fclose(fileu);
	free(fik); fik = NULL;
	free(data); data = NULL;

}


/* Main program */
int main()
{
	//t1();
	//t2();
	t3();
	//t4();
}
