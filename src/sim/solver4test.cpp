#include "solver4.h"
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
    const int m = 1;                    // viscothermal loss IIR filter order
    const float x0 = -2;                 // position of leftmost boundary
    const float x1 = 2;                 // position of rightmost boundary
    const float xmid = (x0 + x1) / 2;   // middle of domain along x
    const float dx=(x1-x0)/(nx-1);      // cell width
    const float dt = 1/4.41e4;          // timestep per iteration
    const float c = 343;                // speed of sound
    const float rho = 1.293;            // density of air
    const float dg = 0.1 * (x1 - x0);   // standard deviation of density array gaussian

    const float S0 = 0.5;
    const float S1 = 1;
    const float dS = (S1-S0)/(nx-1);

    const long sstep = 1e0;

    // initialise arrays
    float x[nx], p[nx], v[nx-1], S[nx];
    auto pdata = new float[(int) (nt / sstep)][nx+2];
    auto vdata = new float[(int) (nt / sstep)][nx+1];

    // fill x array with nx evenly-spaced positions between x0 and x1
    generate(x, x+nx, [i=0, dx, x0]() mutable {return x0 + dx * i++;});

    // fill pressure array with a gaussian blob centered on xmid
    // generate(p, p+nx, [i=0, x, xmid, dg]() mutable {return (float) 1e4 * exp(-pow((x[i++]-xmid)/dg,2));});
    fill(p, p+nx, 0);

    // fill velocity array with const zero
    fill(v, v+nx-1, 0);

    // fill area array with constant 1 area
    
    generate(S, S+nx, [i=0, dS, S0]() mutable {return S0 + dS * i++;});

    // riemann condition

    Solver4<nx> s = Solver4<nx>(dx, dt, c, rho, p, v, S, BoundaryType::Driver, BoundaryType::Open);

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

    ofstream myfile_p;
    myfile_p.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/pdata.bin", ios::binary);
    for (int i = 0; i < nt / sstep; i++) {
        myfile_p.write((char*) pdata[i], sizeof(float)*(nx+2));
    }
    
    myfile_p.close();

    ofstream myfile_v;
    myfile_v.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/vdata.bin", ios::binary);
    for (int i = 0; i < nt / sstep; i++) {
        myfile_v.write((char*) vdata[i], sizeof(float)*(nx+1));
    }
    
    myfile_v.close();

    return 0;
    
}