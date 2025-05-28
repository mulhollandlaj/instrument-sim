#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <iomanip>
#include "solver.h"
#include <fstream>
#include <chrono>

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
    const int nt = 1e5;                // no. iterations
    const float x0 = .3;                 // position of leftmost boundary
    const float x1 = .7;                 // position of rightmost boundary
    const float xmid = (x0 + x1) / 2;   // middle of domain along x
    const float dx=(x1-x0)/(nx-1);      // cell width
    const float dt = 0.005;              // timestep per iteration
    const float cfl = 0.5;              // Courant number
    const float gamma = 1.4;            // ratio of specific heats
    const float dg = 0.1 * (x1 - x0);   // standard deviation of density array gaussian

    // initialise arrays
    float x[nx], rho[nx], rhou[nx], e[nx], time[nt+1];
    auto pdata = new float[nt][nx];

    // fill x array with nx evenly-spaced positions between x0 and x1
    generate(x, x+nx, [i=0, dx, x0]() mutable {return x0 + dx * i++;});

    // fill density array with a gaussian blob centered on xmid
    generate(rho, rho+nx, [i=0, x, xmid, dg]() mutable {return (float) (1 + 0.3*exp(-pow((x[i++]-xmid)/dg,2)));});

    // fill momentum density array with constant 0
    fill(rhou, rhou+nx, 0);

    // fill internal energy array with constant 1
    fill(e, e+nx, 1);

    Solver<nx> s = Solver<nx>(x, rho, rhou, e, gamma);
    auto start = chrono::high_resolution_clock::now();
    for (int i = 1; i < nt+1; i++) {
        s.step(dt);
        s.getPressure(pdata[i-1]);
    }
    auto end = chrono::high_resolution_clock::now();

    double timeTaken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    timeTaken *= 1e-9;

    double timePerIter = timeTaken / (double) nt;

    cout << "Solved " << nt << "iterations in " << timeTaken << setprecision(9) << " s" << endl;

    cout << timePerIter * 1e6 << " us / iteration" << endl;
    cout << 1 / timePerIter << " Hz" << endl;

    // printFloatArr(pdata[0], nx);

    ofstream myfile;
    myfile.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/pdata.bin", ios::binary);
    for (int i = 0; i < nt; i++) {
        myfile.write((char*) pdata[i], sizeof(float)*nx);
    }
    
    myfile.close();

    return 0;
}