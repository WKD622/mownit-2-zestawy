#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

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

double _abs(double x)
{
    if (x < 0.0)
        return x * (-1.0);
    else
        return x;
}

void test(int min_number_of_shoots, int max_number_of_shoots, FILE *results1, FILE *results2)
{
    double x_result = 0.33333333333;
    double sqrt_result = 2.0;
    for (int i = min_number_of_shoots; i < max_number_of_shoots; i++)
    {
        printf("%d\n", i);
        fprintf(results1, "%3d   %10f\n", i, _abs(x_result - hit_and_miss_x_2(i, 0.0, 1.0, 0.0, 1.0)));
        fprintf(results2, "%3d   %10f\n", i, _abs(sqrt_result - hit_and_miss_1_sqr(i, 0.0, 1.0, 0.0, 10000.0)));
    }
}

int main()
{
    srand(time(NULL));
    //FILE *results1 = fopen("out/results_x_2", "wr");
    //FILE *results2 = fopen("out/results_1_sqrt", "wr");

    //test(100, 100000, results1, results2);
    system("gnuplot --persist -e 'plot \"out/results_x_2\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/results_1_sqrt\" u 1:2'");

    //fclose(results1);
    //fclose(results2);

    return 0;
}