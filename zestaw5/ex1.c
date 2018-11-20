#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#define N 10

double x_2(double x)
{
    return pow(x, 2.0);
}

double sqr_1(double x)
{
    return 1.0 / sqrt(x);
}

double rand_value(double start, double end)
{
    return (double)rand() / RAND_MAX * (end - start) + start;
}

double hit_and_miss_x_2(int number_of_shoots, double x_start, double x_end, double y_start, double y_end)
{
    int success = 0;
    int loss = 0;
    for (int i = 0; i < number_of_shoots; i++)
    {
        double x = rand_value(x_start, x_end);
        double y = rand_value(y_start, y_end);

        double f_value = x_2(x);
        if (y <= f_value)
            success++;
        else
            loss++;
    }
    return ((double)success / ((double)loss + (double)success)) * (x_end - x_start) * (y_end - y_start);
}

double hit_and_miss_1_sqr(int number_of_shoots, double x_start, double x_end, double y_start, double y_end)
{
    int success = 0;
    int loss = 0;
    for (int i = 0; i < number_of_shoots; i++)
    {
        double x = rand_value(x_start, x_end);
        double y = rand_value(y_start, y_end);

        double f_value = sqr_1(x);
        if (y <= f_value)
            success++;
        else
            loss++;
    }
    return ((double)success / ((double)loss + (double)success)) * (x_end - x_start) * (y_end - y_start);
}

int main()
{
    srand(time(NULL));
    printf("%f\n", hit_and_miss_x_2(100000, 0.0, 1.0, 0.0, 1.0));
    printf("%f\n", hit_and_miss_1_sqr(100000000, 0.0, 1.0, 0.0, 10000.0));
    return 0;
}