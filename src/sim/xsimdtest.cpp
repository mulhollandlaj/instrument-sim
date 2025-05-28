#include "xsimd/xsimd.hpp"
#include <stdio.h>
#include <iostream>


namespace xs = xsimd;
const int BSIZE = sizeof(xs::batch<float>) / sizeof(float);


xs::batch<float> mean(xs::batch<float> lhs, xs::batch<float> rhs) {
    return (lhs + rhs) / 2;
}

int main() {
    std::cout << BSIZE << std::endl;

    // defining constants for problem
    const int nx = 512;
    const int nt = 1e4;
    const float x0 = 0;
    const float x1 = 1;
    const float xmid = (x0 + x1) / 2;
    const float dx = (x1-x0)/(nx-1);
    const float dt = 0.0025;
    const float cfl = 0.5;
    const float gamma = 1.4;
    const float dg = 0.1 * (x1 - x0);

    const int nb = nx / BSIZE;

    xs::batch<float> xl[nb], x[nb], xr[nb], rhol[nb], rho[nb], rhor[nb], rhoul[nb], rhou[nb], rhour[nb];
}
