#include "solver2.h"
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace std;

int main() {
    // defining constants for problem
    const int nx = 128;                 // no. cells along x
    const long nt = 1e5;                // no. iterations
    const float x0 = 0;                 // position of leftmost boundary
    const float x1 = 3;                 // position of rightmost boundary
    const float xmid = (x0 + x1) / 2;   // middle of domain along x
    const float dx=(x1-x0)/(nx-1);      // cell width
    const float dt = 5e-6;             // timestep per iteration
    const float cfl = 0.5;              // Courant number
    const float gamma = 1.4;            // ratio of specific heats
    const float dg = 0.1 * (x1 - x0);   // standard deviation of density array gaussian

    const float S0 = 0.1;
    const float S1 = 1;
    const float dS = (S1-S0)/(nx-1);

    const long sstep = 1e1;

    const int ampl = 1;
    const int freq = 100;

    // initialise arrays
    float x[nx], S[nx], rho[nx], rhou[nx], e[nx];
    auto pdata = new float[(int) (nt / sstep)][nx+2];
    auto vdata = new float[(int) (nt / sstep)][nx+1];

    // fill x array with nx evenly-spaced positions between x0 and x1
    generate(x, x+nx, [i=0, dx, x0]() mutable {return x0 + dx * i++;});

    // fill area array with constant 1 area
    // fill(rhou, rhou+nx, 1);
    generate(S, S+nx, [i=0, dS, S0]() mutable {return S0 + dS * i++;});

    // fill density array with a gaussian blob centered on xmid
    // generate(rho, rho+nx, [i=0, x, xmid, dg, S]() mutable {return (float) 1.293 * (1 + 0.1*exp(-pow((x[i]-xmid)/dg,2))) /S[i++];});
    fill(rho, rho+nx/2,1.293);
    fill(rho+nx/2, rho+nx, 1.293*1.1);

    // riemann condition
    // generate(rho, rho+nx/2, [i=0, S]() mutable {return (float) 1.293 * 1.2 / S[i++];});
    // generate(rho+nx/2, rho+nx, [i=nx/2, S]() mutable {return (float) 1.293 / S[i++];});

    // fill momentum density array with constant 0
    fill(rhou, rhou+nx, 0);

    // fill internal energy array with constant 1
    fill(e, e+nx, 214e3);

    Solver3<nx> s = Solver3<nx>(x, S, rho, rhou, e, dt, gamma);
    s.setBoundaries(rho[0], rho[nx-1], rhou[0], rhou[nx-1]);

    auto start = chrono::high_resolution_clock::now();
    for (long i = 0; i < nt / sstep; i++) {
        for (long j = 0; j < sstep; j++) {
            // s.setBoundaries(1 + ampl * sin(2*M_PI*freq*s.getTime()), rho[nx-1], 0, 0);
            s.step();
        }
        s.getPressure(pdata[i]);
        s.getVelocity(vdata[i]);
    }
    auto end = chrono::high_resolution_clock::now();

    double timeTaken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    timeTaken *= 1e-9;

    double timePerIter = timeTaken / (double) nt;

    cout << "Solved " << nt << " iterations in " << timeTaken << setprecision(9) << " s" << endl;

    cout << timePerIter * 1e6 << " us / iteration" << endl;
    cout << 1 / timePerIter << " Hz" << endl;

    ofstream myfilep;
    myfilep.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/pdata.bin", ios::binary);
    for (int i = 0; i < nt / sstep; i++) {
        myfilep.write((char*) pdata[i], sizeof(float)*(nx+2));
    }

    myfilep.close();

    ofstream myfilev;
    myfilev.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/vdata.bin", ios::binary);
    for (int i = 0; i < nt / sstep; i++) {
        myfilev.write((char*) vdata[i], sizeof(float)*(nx+1));
    }
    
    myfilev.close();

    return 0;
    
}