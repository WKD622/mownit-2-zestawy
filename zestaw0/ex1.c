#include <gsl/gsl_interp.h> 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#define N 100

float float_1[N];
double double_1[N];
float float_2[N];
double double_2[N];
double diff_1[N];
double diff_2[N];

void prepare()
{
    int i;
    for (i = 0; i < N; i++){
        float_1[i] = 0;
        double_1[i] = 0;
        float_2[i] = 0;
        double_2[i] = 0;
    }
}

float x_float_1(int n)
{   
    float result;
    if (n==0) return 0.01;
    if (float_1[n] != 0) 
        return float_1[n];
    result = x_float_1(n-1) + 3.0 * x_float_1(n-1) * (1-x_float_1(n-1));
    if (float_1[n] == 0) 
        float_1[n] = result;
    return result;
}

double x_double_1(int n)
{
    double result;
    if (n==0) return 0.01;
    if (double_1[n] != 0) 
        return double_1[n];
    result = x_double_1(n-1) + 3.0 * x_double_1(n-1) * (1-x_double_1(n-1));
    if (double_1[n] == 0)
        double_1[n] = result;
    return result;
}

float x_float_2(int n)
{
    float result;
    if (n==0) return 0.01;
    if (float_2[n] != 0) 
        return float_2[n];
    result =  4.0 * x_float_2(n-1) - 3.0 * x_float_2(n-1) * x_float_2(n-1);
    if (float_2[n] == 0) 
        float_2[n] = result;
    return result;
}

double x_double_2(int n)
{
    double result;
    if (n==0) return 0.01;
    if (double_2[n] != 0) 
        return double_2[n];
    result =  4.0 * x_double_2(n-1) - 3.0 * x_double_2(n-1) * x_double_2(n-1);
    if (double_2[n] == 0)
        double_2[n] = result;
    return result;
}

double mod(double l)
{
    if (l < 0)
        return l*(-1);
    else return l;
}

void save_to_files()
{
    int i;
    FILE *diff_1_f = fopen("out/diff_1_f","w+");
    FILE *diff_2_f = fopen("out/diff_2_f","w+");
    FILE *x_float_1_f = fopen("out/x_float_1_f","w+");
    FILE *x_double_1_f = fopen("out/x_double_1_f","w+");
    FILE *x_float_2_f = fopen("out/x_float_2_f","w+");
    FILE *x_double_2_f = fopen("out/x_double_2_f","w+");
    for (i = 0; i < N; i++){
        fprintf(diff_1_f, "%i,%f\n", i, diff_1[i]);
        fprintf(diff_2_f, "%i,%f\n", i, diff_2[i]);
        fprintf(x_float_1_f, "%i,%f\n", i, float_1[i]);
        fprintf(x_double_1_f, "%i,%f\n", i, double_1[i]);
        fprintf(x_float_2_f, "%i,%f\n", i, float_2[i]);
        fprintf(x_double_2_f, "%i,%f\n", i, double_2[i]);
    }
}

double epsilon(){
    double e = 1;
    while (e + 1.0 > 1.0)
        e = e/2.0;
    return e;
}

int main()
{
    prepare();
    int i;
    for (i = 0; i < N; i++){
        float x_float_result = x_float_1(i);
        double x_double_result = x_double_1(i);
        double diff = mod(x_double_result - (double)x_float_result);
        diff_1[i] = diff;
        //printf("[1st] for x(%d):\nfloat:  %.20f\ndouble: %.20f\ndiffrence: %.20f\n\n", i, x_float_result, x_double_result, diff);
    }

    for (i = 0; i < N; i++){
        float x_float_result = x_float_2(i);
        double x_double_result = x_double_2(i);
        double diff = mod(x_double_result - (double)x_float_result);
        diff_2[i] = diff;
        //printf("[2nd] for x(%d):\nfloat:  %.20f\ndouble: %.20f\ndiffrence: %.20f\n\n", i, x_float_result, x_double_result, diff);
    }
    save_to_files();

    printf("%.30f", epsilon());
    return 0;
}
