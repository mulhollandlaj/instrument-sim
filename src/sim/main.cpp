#include <immintrin.h>
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "solver.h"
#include <fstream>

using namespace std;

void printFloatArr(float* arr, int size) {
    for (int i = 0; i < size; i++) {
        printf( "%6.*lf ", 4, arr[i]);
    }
    cout << endl;
}

int main() {
    // defining constants for problem
    const int nx = 128;                 // no. cells along x
    const int nt = 1000;                // no. iterations
    const float x0 = 0;                 // position of leftmost boundary
    const float x1 = 1;                 // position of rightmost boundary
    const float xmid = (x0 + x1) / 2;   // middle of domain along x
    const float dx=(x1-x0)/(nx-1);      // cell width
    const float dt = 0.0025;              // timestep per iteration
    const float cfl = 0.5;              // Courant number
    const float gamma = 1.4;            // ratio of specific heats
    const float dg = 0.1 * (x1 - x0);   // standard deviation of density array gaussian

    // initialise arrays
    float x[nx], rho[nx], rhou[nx], e[nx], time[nt+1], pdata[nt][nx];

    // fill x array with nx evenly-spaced positions between x0 and x1
    generate(x, x+nx, [i=0, dx, x0]() mutable {return x0 + dx * i++;});

    // fill density array with a gaussian blob centered on xmid
    generate(rho, rho+nx, [i=0, x, xmid, dg]() mutable {return (float) (1 + 0.3*exp(-pow((x[i++]-xmid)/dg,2)));});

    // fill momentum density array with constant 0
    fill(rhou, rhou+nx, 0);

    // fill internal energy array with constant 1
    fill(e, e+nx, 1);

    Solver<128> s = Solver<128>(x, rho, rhou, e, gamma);
    for (int i = 1; i < nt+1; i++) {
        s.hydroIso(dt);
        s.getPressure(pdata[i-1]);
    }

    printFloatArr(pdata[0], nx);

    cout << "Hello, world!" << endl;

    ofstream myfile;
    myfile.open("pdata.bin", ios::binary);
    myfile.write((char*) &pdata, sizeof(pdata));
    myfile.close();

    return 0;
}


void avxDemo() {
    /* Initialize the two argument vectors */
    __m256 evens = _mm256_set_ps(2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0);
    __m256 odds = _mm256_set_ps(1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0);

    /* Compute the difference between the two vectors */
    __m256 result = _mm256_sub_ps(evens, odds);

    /* Display the elements of the result vector */
    float *f = (float *)&result;
    printf("%f %f %f %f %f %f %f %f\n",
          f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
}