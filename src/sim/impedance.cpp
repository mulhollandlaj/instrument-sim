#include "solver4.h"
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace std;

int main() {
    const int nx = 128;                 // no. cells along x
    const long nt = 1e5;                // no. iterations
    const int m = 13;                   // viscothermal loss fractional derivative approximation order
    const float x0 = -2;                 // position of leftmost boundary
    const float x1 = 2;                 // position of rightmost boundary
    const float xmid = (x0 + x1) / 2;   // middle of domain along x
    const float dx=(x1-x0)/(nx-1);      // cell width
    const float dt = 1/4.41e4;          // timestep per iteration
    const float c = 3.4723e2;           // speed of sound
    const float rho = 1.1769;           // density of air
    const float eta = 1.846e-5;         // shear viscosity coefficient
    const float nu = 0.8410;            // root Prandtl number
    const float gamma = 1.4017;         // ratio of specific heats

    const float dg = 0.1 * (x1 - x0);   // standard deviation of density array gaussian

    const float a0 = 0.05;
    const float a1 = 0.1;
    const float da = (a1-a0)/(nx-1);

    const long sstep = 1e2;

    // initialise arrays
    float x[nx], p[nx+1], v[nx-1], a[nx], S[nx];
    float pdata[nt];
    float vdata[nt];

    // fill x array with nx evenly-spaced positions between x0 and x1
    generate(x, x+nx, [i=0, dx, x0]() mutable {return x0 + dx * i++;});

    // fill pressure array with a gaussian blob centered on xmid
    // generate(p, p+nx, [i=0, x, xmid, dg]() mutable {return (float) 1e4 * exp(-pow((x[i++]-xmid)/dg,2));});

    // pressure riemann condition
    // fill(p, p+nx/2, 1e3);
    // fill(p+nx/2, p+nx, -1e3);

    // pressure flat zero
    fill(p, p+nx, 0);

    // fill velocity array with const zero
    fill(v, v+nx-1, 0);

    // fill area array with constant 1 area
    // fill(rhou, rhou+nx, 1);
    generate(a, a+nx, [i=0, da, a0]() mutable {return a0 + da * i++;});

    generate(S, S+nx, [i=0, a]() mutable {return M_PI * pow(a[i++],2);});

    Solver4<nx, m> s = Solver4<nx, m>(dx, dt, c, rho, p, v, S, BoundaryType::Driver, BoundaryType::Closed);

    auto start = chrono::high_resolution_clock::now();
    for (long i = 0; i < nt / sstep; i++) {
        for (long j = 0; j < sstep; j++) {
            // s.setBoundaries(1 + ampl * sin(2*M_PI*freq*s.getTime()), rho[nx-1], 0, 0);
            s.step();
        }
        pdata[i] = s.getPressure(0);
        vdata[i] = s.getVelocity(0);
    }
    auto end = chrono::high_resolution_clock::now();

    double timeTaken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    timeTaken *= 1e-9;

    double timePerIter = timeTaken / (double) nt;

    cout << "Solved " << nt << " iterations in " << timeTaken << setprecision(9) << " s" << endl;

    cout << timePerIter * 1e6 << " us / iteration" << endl;
    cout << 1 / timePerIter << " Hz" << endl;

    float p_rms, v_rms = 0;

    for(int i = (int) (0.1 / dt); i < nt; i++) {
        p_rms += pow(pdata[i],2);
        v_rms += pow(vdata[i],2);
    }
    p_rms = sqrt(p_rms) / nt;
    v_rms = sqrt(v_rms) / nt;

    cout << "RMS Pressure: " << p_rms << endl << "RMS Velocity: " << v_rms << endl << "Impedance: " << p_rms / v_rms << endl;
}