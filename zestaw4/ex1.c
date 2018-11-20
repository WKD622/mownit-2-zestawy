#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#define N 256
#define Pi 3.14159265358979323846

void print_array(double T[N])
{
    int i;
    for (i = 0; i < N; i++)
    { 
        printf("%d) %10f\n", i, T[i]);
    }
}

void generate_array1(double *T)
{
    int i;
    for (i = 0; i < N; i++)
    {
        T[i] = cos(4.0 * Pi * (double)i / (double)N) + cos(16.0 * Pi * (double)i / (double)N) / 5.0 + cos(32.0 * Pi * (double)i / (double)N) / 8.0 + cos(128.0 * Pi * (double)i / (double)N) / 16.0;
    }
}

void generate_array2(double *T)
{
    int i;
    for (i = 0; i < N; i++)
    {
        T[i] = cos(4.0 * Pi * (double)i / (double)N) + ((float)rand()) / RAND_MAX / 8.0;
    }
}

double _abs(double l)
{
    if (l < 0.0)
        return l * (-1.0);
    else
        return l;
}

void generate_file(FILE *results, double *T)
{
    int i;
    for (i = 0; i < N; i++)
    {
        fprintf(results, "%3d   %10f\n", i, T[i]);
    }
}

void generate_file_abs(FILE *results, double *T)
{
    int i;
    for (i = 0; i < N; i++)
    {
        if (_abs(T[i]) < 50)
            fprintf(results, "%3d   %10f\n", i, 0.0);
        else
            fprintf(results, "%3d   %10f\n", i, T[i]);
    }
}

void clear_function(double *T)
{
    int i;
    for (i = 0; i < N; i++)
    {
        if (_abs(T[i]) < 10)
            T[i] = 0.0;
    }
}

int main()
{
    double *T = malloc(N * sizeof(double));
    double *S = malloc(N * sizeof(double));
    srand(time(NULL));
    generate_array1(T);
    generate_array2(S);
    FILE *results1 = fopen("out/results1", "wr");
    FILE *results2 = fopen("out/results2", "wr");
    FILE *results3 = fopen("out/results3", "wr");
    FILE *results4 = fopen("out/results4", "wr");
    FILE *results5 = fopen("out/results5", "wr");
    FILE *results6 = fopen("out/results6", "wr");
    generate_file(results1, T);
    generate_file(results2, S);

    gsl_fft_real_radix2_transform(T, 1, 256);
    gsl_fft_real_radix2_transform(S, 1, 256);
    generate_file(results3, T);
    generate_file_abs(results4, S);

    clear_function(S);
    clear_function(T);
    gsl_fft_halfcomplex_radix2_inverse(S, 1, 256);
    gsl_fft_halfcomplex_radix2_inverse(T, 1, 256);
    generate_file(results5, S);
    generate_file(results6, T);

    system("gnuplot --persist -e 'plot \"out/results1\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/results3\" u 1:2 with boxes'");
    system("gnuplot --persist -e 'plot \"out/results2\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/results4\" u 1:2 with boxes'");
    system("gnuplot --persist -e 'plot \"out/results5\" u 1:2'");
    //system("gnuplot --persist -e 'plot \"out/results6\" u 1:2'");
    free(T);
    free(S);
    fclose(results1);
    fclose(results2);
    fclose(results3);
    fclose(results4);
    fclose(results5);
    fclose(results6);
}
