#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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

int main()
{
    double T[N];
    int i;
    for (i = 0; i < N; i++)
    {
        T[i] = cos(4.0 * Pi * (double)i / (double)N) + cos(16.0 * Pi * (double)i / (double)N) / 5.0 + cos(32.0 * Pi * (double)i / (double)N) / 8.0 + cos(128.0 * Pi * (double)i / (double)N) / 16.0;
    }
    print_array(T);
}