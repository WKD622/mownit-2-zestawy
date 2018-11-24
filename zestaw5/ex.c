#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

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

void test_my_monte_carlo(int min_number_of_shoots, int max_number_of_shoots, FILE *results1, FILE *results2)
{
    double x_result = 0.3333333333333;
    double sqrt_result = 2.0;
    for (int i = min_number_of_shoots; i < max_number_of_shoots; i++)
    {
        printf("my monte carlo: %d of %d\n", i, max_number_of_shoots);
        system("clear");
        fprintf(results1, "%3d   %30.25f\n", i, _abs(x_result - hit_and_miss_x_2(i, 0.0, 1.0, 0.0, 1.0)));
        fprintf(results2, "%3d   %30.25f\n", i, _abs(sqrt_result - hit_and_miss_1_sqr(i, 0.0, 1.0, 0.0, 10000.0)));
    }
}

double x_2_gsl(double *x, size_t dim, void *p)
{
    return x[0] * x[0];
}

void test_x_2_gsl(int min_number_of_shoots, int max_number_of_shoots, FILE *plain, FILE *miser, FILE *vegas)
{
    double x_result = 0.3333333333333;
    for (int i = min_number_of_shoots; i < max_number_of_shoots; i++)
    {
        double res, err;
        double xl[2] = {0};
        double xu[2] = {1};
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_monte_function G = {&x_2_gsl, 1, 0};
        size_t calls = i;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_monte_plain_state *s1 = gsl_monte_plain_alloc(1);
        gsl_monte_plain_integrate(&G, xl, xu, 1, calls, r, s1,
                                  &res, &err);
        gsl_monte_plain_free(s1);
        fprintf(plain, "%3d   %30.25f\n", i, _abs(x_result - res));

        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_monte_miser_state *s2 = gsl_monte_miser_alloc(1);
        gsl_monte_miser_integrate(&G, xl, xu, 1, calls, r, s2,
                                  &res, &err);
        gsl_monte_miser_free(s2);
        fprintf(miser, "%3d   %30.25f\n", i, _abs(x_result - res));

        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_monte_vegas_state *s3 = gsl_monte_vegas_alloc(1);
        gsl_monte_vegas_integrate(&G, xl, xu, 1, calls, r, s3,
                                  &res, &err);
        fprintf(vegas, "%3d   %30.25f\n", i, _abs(x_result - res));
        gsl_monte_vegas_free(s3);

        printf("x^2 gsl: %d of %d\n", i, max_number_of_shoots);
        system("clear");
    }
}

double sqrt_1_gsl(double *x, size_t dim, void *p)
{
    return 1.0 / sqrt(x[0]);
}

void test_sqrt_1_gsl(int min_number_of_shoots, int max_number_of_shoots, FILE *plain, FILE *miser, FILE *vegas)
{
    double sqrt_result = 2.0;
    for (int i = min_number_of_shoots; i < max_number_of_shoots; i++)
    {
        double res, err;
        double xl[2] = {0};
        double xu[2] = {1};
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_monte_function G = {&sqrt_1_gsl, 1, 0};
        size_t calls = i;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_monte_plain_state *s1 = gsl_monte_plain_alloc(1);
        gsl_monte_plain_integrate(&G, xl, xu, 1, calls, r, s1,
                                  &res, &err);
        gsl_monte_plain_free(s1);
        fprintf(plain, "%3d   %30.25f\n", i, _abs(sqrt_result - res));

        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_monte_miser_state *s2 = gsl_monte_miser_alloc(1);
        gsl_monte_miser_integrate(&G, xl, xu, 1, calls, r, s2,
                                  &res, &err);
        gsl_monte_miser_free(s2);
        fprintf(miser, "%3d   %30.25f\n", i, _abs(sqrt_result - res));

        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_monte_vegas_state *s3 = gsl_monte_vegas_alloc(1);
        gsl_monte_vegas_integrate(&G, xl, xu, 1, calls, r, s3,
                                  &res, &err);
        fprintf(vegas, "%3d   %30.25f\n", i, _abs(sqrt_result - res));
        gsl_monte_vegas_free(s3);

        printf("1/sqrt(x) gsl: %d of %d\n", i, max_number_of_shoots);
        system("clear");
    }
}

int main()
{
    srand(time(NULL));

    // testing my monte carlo for x^2 and 1/sqrt(x)
    FILE *my_x_2 = fopen("out/my_x_2", "wr");
    FILE *my_1_sqrt = fopen("out/my_1_sqrt", "wr");
    test_my_monte_carlo(100, 50000, my_x_2, my_1_sqrt);
    fclose(my_x_2);
    fclose(my_1_sqrt);

    // testing gsl monte carlo for x^2
    FILE *plain_x_2 = fopen("out/plain_x_2", "wr");
    FILE *miser_x_2 = fopen("out/miser_x_2", "wr");
    FILE *vegas_x_2 = fopen("out/vegas_x_2", "wr");
    test_x_2_gsl(100, 50000, plain_x_2, miser_x_2, vegas_x_2);
    fclose(plain_x_2);
    fclose(miser_x_2);
    fclose(vegas_x_2);

    // testing gsl monte carlo for 1/sqrt(x)
    FILE *plain_sqrt_1 = fopen("out/plain_sqrt_1", "wr");
    FILE *miser_sqrt_1 = fopen("out/miser_sqrt_1", "wr");
    FILE *vegas_sqrt_1 = fopen("out/vegas_sqrt_1", "wr");
    test_sqrt_1_gsl(100, 50000, plain_x_2, miser_x_2, vegas_x_2);
    fclose(plain_sqrt_1);
    fclose(miser_sqrt_1);
    fclose(vegas_sqrt_1); 

    // drawing
    system("gnuplot --persist -e 'plot \"out/my_x_2\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/my_1_sqrt\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/plain_x_2\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/miser_x_2\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/vegas_x_2\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/plain_sqrt_1\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/miser_sqrt_1\" u 1:2'");
    system("gnuplot --persist -e 'plot \"out/vegas_sqrt_1\" u 1:2'");

    return 0;
}
