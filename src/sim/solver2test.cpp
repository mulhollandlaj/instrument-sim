#include "solver2.h"
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace std;

int main() {
    // defining constants for problem
    const int nx = 512;                 // no. cells along x
    const long nt = 1e5;                // no. iterations
    const float x0 = 0;                 // position of leftmost boundary
    const float x1 = 1;                 // position of rightmost boundary
    const float xmid = (x0 + x1) / 2;   // middle of domain along x
    const float dx=(x1-x0)/(nx-1);      // cell width
    const float dt = 1./4.4e4;              // timestep per iteration
    const float cfl = 0.5;              // Courant number
    const float gamma = 1.4;            // ratio of specific heats
    const float dg = 0.1 * (x1 - x0);   // standard deviation of density array gaussian

    const float A0 = 0.2;
    const float A1 = 1.2;
    const float dA = (A1-A0)/(nx-1);

    const long sstep = 1e3;

    const int ampl = 1;
    const int freq = 100;

    // initialise arrays
    float x[nx], A[nx], rho[nx], rhou[nx], rhoe[nx];
    auto pdata = new float[(int) (nt / sstep)][nx];

    // fill x array with nx evenly-spaced positions between x0 and x1
    generate(x, x+nx, [i=0, dx, x0]() mutable {return x0 + dx * i++;});

    // fill area array with constant 1 area
    generate(A, A+nx, [i=0, dA, A0]() mutable {return A0 + dA * i++;});

    // fill density array with a gaussian blob centered on xmid
    generate(rho, rho+nx, [i=0, x, xmid, dg]() mutable {return (float) 10 * (1 + 0.3*exp(-pow((x[i++]-xmid)/dg,2)));});

    // fill momentum density array with constant 0
    fill(rhou, rhou+nx, 0);

    // fill internal energy array with constant 1
    // copy(rho, rho+nx, rhoe);
    generate(rhoe, rhoe+nx, [i=0, x, xmid, dg]() mutable {return (float) 1 * (1 + 0.3*exp(-pow((x[i++]-xmid)/dg,2)));});

    Solver2<nx> s = Solver2<nx>(x, A, rho, rhou, rhoe, dt, gamma);
    s.setBoundaries(rho[0], rho[nx-1], 0, 0, rhoe[0], rhoe[nx-1]);

    auto start = chrono::high_resolution_clock::now();
    for (long i = 0; i < nt / sstep; i++) {
        for (long j = 0; j < sstep; j++) {
            // s.setBoundaries(1 + ampl * sin(2*M_PI*freq*s.getTime()), rho[nx-1], 0, 0, rhoe[0], rhoe[nx-1]);
            s.step();
        }
        s.getPressure(pdata[i]);
    }
    auto end = chrono::high_resolution_clock::now();

    double timeTaken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    timeTaken *= 1e-9;

    double timePerIter = timeTaken / (double) nt;

    cout << "Solved " << nt << "iterations in " << timeTaken << setprecision(9) << " s" << endl;

    cout << timePerIter * 1e6 << " us / iteration" << endl;
    cout << 1 / timePerIter << " Hz" << endl;

    ofstream myfile;
    myfile.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/pdata.bin", ios::binary);
    for (int i = 0; i < nt / sstep; i++) {
        myfile.write((char*) pdata[i], sizeof(float)*nx);
    }
    
    myfile.close();

    return 0;
    
}