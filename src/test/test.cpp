#define FULL_OUT
#define WAV_OUT
// #define DEBUG

#include <algorithm>
#include <fstream>
#include <chrono>
#include <iostream>
#include <iomanip>

#include "../sim2/solver.h"
#include "../sim2/webster.h"
#include "../sim2/visco.h"
#include "../sim2/euler.h"
#include "../test/wavwriter.h"
#include "../test/profile.h"

template <int nx>
void genHorn(std::array<float, nx> &SArr, float l, float a_0, float a_l, int n)
{
    float h = l / (nx - 1);
    float c_0 = std::pow(a_0, 1.0 / n);
    float c_l = std::pow(a_l, 1.0 / n);
    std::generate(SArr.begin(), SArr.end(), [l, h, c_0, c_l, n, i = 0]() mutable
                  { return std::pow((c_l - c_0) / l * h * i++ + c_0, n); });
}

template <int nx>
void genRiemann(std::array<float, nx> &arr, float nval, float pval) {
    std::fill(arr.begin(), arr.begin()+nx/2, nval);
    std::fill(arr.begin()+nx/2, arr.end(), pval);
}

int main()
{
    const int nx = 64;
    const int nt = 8.82e4;
    const int m = 12;
    const float sim_freq = 8.82e4;

    profile::Profile pr = profile::Profile<nx>("C:\\Users\\mulhollandlaj\\Documents\\4YP\\instrument-sim\\input\\delusse_oboe.bin");

    const float driver_ampl = 1;
    const float driver_freq = 110;

    const int bufferSize = 256;

    float a_0 = 2e-3;
    float a_l = 2e-1;
    int n = 6;

    const float k = 1.0 / sim_freq;
    const float h = pr.length / (nx - 1);
    const float c = 3.4723e2;   // speed of sound
    const float rho = 1.1769;   // density of air
    const float eta = 1.846e-5; // shear viscosity coefficient
    const float nu = 0.8410;    // root Prandtl number
    const float gamma = 1.4017; // ratio of specific heats
    const float e0 = 214e3;
    const float p0 = 1.293 * (1.4 - 1) * e0;

    std::array<float, nx> p, rho_, rho_u_, e;
    std::array<float, nx + 1> v;

    std::fill(p.begin(), p.end(), 0);
    std::fill(v.begin(), v.end(), 0);
    v[1] = 1 / h;
    // p[0] = h * rho / k * v[0];

    // std::generate(S.begin(), S.end(), [S0, SL, i=0]() mutable {return S0 + i++ * (SL - S0) / (nx-1);});
    // genHorn<S.size()>(S, l, a_0, a_l, n);

    std::generate(rho_.begin(), rho_.end(), [p, pr, c, i = 0]() mutable {return (p[i]*std::pow(c,-2) / pr.S[i++]);});

    // fill momentum density array with constant 0
    std::fill(rho_u_.begin(), rho_u_.end(), 0);

    // fill internal energy array with constant 1
    std::fill(e.begin(), e.end(), e0);

    std::array<float, nx> pDest, vDest;

    std::array<short, bufferSize> pBuffer;

    using WS = solver::WebsterSolver<nx>;
    WS ws = WS(
        k,
        h,
        rho,
        c,
        p,
        v,
        pr.S,
        solver::BFuncClosed<WS>(),
        solver::BFuncClosed<WS>()
    );

    using VS = solver::ViscoSolver<nx, m>;
    VS vs = VS(
        k,
        h,
        rho,
        c,
        eta,
        nu,
        gamma,
        p,
        v,
        pr.S,
        solver::BFuncClosed<VS>(),
        solver::BFuncClosed<VS>()
    );

    using ES = solver::EulerSolver<nx>;
    ES es = ES(
        k,
        h,
        1.4,
        p0,
        rho,
        rho_,
        rho_u_,
        pr.S,
        e,
        solver::BFuncClosed<ES>(),
        solver::BFuncClosed<ES>()
        
    );

    // solver::BFuncDriver<ES>(driver_ampl, driver_freq, 0)
    
    std::cout << "CFL: " << c * k / h << std::endl;

#ifdef WAV_OUT
    WavWriter ww = WavWriter((std::string) "C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/test.wav", 1, 1, (int)sim_freq, 16);
#endif

#ifdef FULL_OUT
    std::ofstream myfile_p;
    std::ofstream myfile_v;
    std::ofstream myfile_s;

    myfile_p.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/pdata.bin", std::ios::binary);
    myfile_v.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/vdata.bin", std::ios::binary);
    myfile_s.open("C:/Users/mulhollandlaj/Documents/4YP/instrument-sim/output/sdata.bin", std::ios::binary);

    
    myfile_p.write((const char *)&nt, sizeof(nt));
    myfile_p.write((const char *)&nx, sizeof(nx));
    myfile_v.write((const char *)&nt, sizeof(nt));
    myfile_v.write((const char *)&nx, sizeof(nx));
    myfile_s.write((const char *)&nx, sizeof(nx));
    myfile_s.write((const char*)&pr.length, sizeof(pr.length));
    myfile_s.write((const char *)&pr.S, nx * sizeof(float));
#endif

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < nt / bufferSize; i++)
    {
        for (int j = 0; j < bufferSize; j++)
        {
            ws.step();
            vs.step();
            es.step();
#ifdef FULL_OUT
            vs.getPressure(pDest);
            vs.getVelocity(vDest);
            
            myfile_p.write((const char *)pDest.data(), nx * sizeof(float));
            myfile_v.write((const char *)vDest.data(), nx * sizeof(float));
#endif

#ifdef WAV_OUT
            float pSample = 1e2 * vs.getPressure(pr.length);
            pBuffer[j] = (short)std::clamp(pSample, (float)-32768, (float)32767);
#endif
        }
#ifdef WAV_OUT
        ww.writeData((char *)pBuffer.data(), bufferSize * sizeof(pBuffer[0]) / sizeof(char));
#endif
    }
    auto end = std::chrono::high_resolution_clock::now();

#ifdef WAV_OUT
    ww.close();
#endif

    double timeTaken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * 1e-9;
    double timePerIter = timeTaken / (double)nt;

    std::cout << "Solved " << nt << " iterations in " << timeTaken << std::setprecision(9) << " s" << std::endl;

    std::cout << timePerIter * 1e6 << " us / iteration" << std::endl;
    std::cout << 1 / timePerIter << " Hz" << std::endl;

#ifdef FULL_OUT
    myfile_p.close();
    myfile_v.close();
#endif
}