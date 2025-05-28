#include <math.h>
#include "xsimd/xsimd.hpp"


namespace xs = xsimd;
const int BSIZE = sizeof(xs::batch<float>) / sizeof(float);


template <int nx>
class Solver2 {
    private:
        float x[nx + 2];
        float dt_dx[nx];

        float A[nx + 2];
        float Anorm[nx];

        float q1[2*(nx+2)];
        float q2[2*(nx+2)];
        float q3[2*(nx+2)];

        float q1i[nx+1];
        float q2i[nx+1];
        float q3i[nx+1];

        float dt, t, c1, c2, c3, c4;
        float *q1o, *q2o, *q3o, *q1n, *q2n, *q3n;

        void swapArrs() {
            bool aOld = q1 == q1o;
            q1o = q1 + (nx+2)*aOld;
            q2o = q2 + (nx+2)*aOld;
            q3o = q3 + (nx+2)*aOld;
            q1n = q1 + (nx+2)*!aOld;
            q2n = q2 + (nx+2)*!aOld;
            q3n = q3 + (nx+2)*!aOld;
        }

    public:
        Solver2(float *x, float *A, float *rho, float *rhou, float *rhoe, float dt, float gamma) : dt(dt) {
            q1o = q1; q2o = q2; q3o = q3;
            q1n = q1+nx+2; q2n = q2+nx+2; q3n = q3+nx+2;

            std::copy(x, x+nx, this->x+1);
            std::copy(A, A+nx, this->A+1);
            std::copy(rho, rho+nx, this->q1o+1);
            std::copy(rhou, rhou+nx, this->q2o+1);
            std::copy(rhoe, rhoe+nx, this->q3o+1);
            
            t = 0;
            c1 = gamma - 1;
            c2 = (3 - gamma) / 2;
            c3 = gamma;
            c4 = (1 - gamma) / 2;

            this->x[0] = 2 * this->x[1] - this->x[2];
            this->x[nx+1] = 2 * this->x[nx] - this->x[nx-1];

            this->A[0] = this->A[1];
            this->A[nx+1] = this->A[nx];

            for (int i = 0; i < nx; i++) {
                dt_dx[i] = 2 * this->dt / (this->x[i+2] - this->x[i]);
                Anorm[i] = (this->A[i+2] - this->A[i]) / (2 * this->A[i+1]);
            }

        }

        void setBoundaries(float q1l, float q1r, float q2l, float q2r, float q3l, float q3r) {
            q1o[0] = q1l; q1n[0] = q1l;
            q1o[nx+1] = q1r; q1n[nx+1] = q1r;
            q2o[0] = q2l; q2n[0] = q2l;
            q2o[nx+1] = q2r; q2n[nx+1] = q2r;
            q3o[0] = q3l; q3n[0] = q3l;
            q3o[nx+1] = q3r; q3n[nx+1] = q3r;
        }

        void step() {
            for (int i = 0; i < nx+1; i++) {                
                q1i[i] = (q1o[i] + q1o[i+1]) / 2;
                q2i[i] = (q2o[i] + q2o[i+1]) / 2;
                q3i[i] = (q3o[i] + q3o[i+1]) / 2;
            }

            for (int i = 0; i < nx; i++) {
                q1n[i+1] = q1o[i+1] - dt_dx[i] * (q2i[i+1] - q2i[i] + Anorm[i] * q2o[i+1]);
                q2n[i+1] = q2o[i+1] - dt_dx[i] * (c1 * (q3i[i+1] - q3i[i]) + c2 * (std::pow(q2i[i+1],2)/q1i[i+1] - std::pow(q2i[i],2)/q1i[i]) + Anorm[i] * (c1 * q3o[i+1] + c2 * std::pow(q2o[i+1],2)/q1o[i+1]));
                q3n[i+1] = q3o[i+1] - dt_dx[i] * (c3 * (q3i[i+1]*q2i[i+1]/q1i[i+1] - q3i[i]*q2i[i]/q1i[i]) + c4 * (std::pow(q2i[i+1],3) / std::pow(q1i[i+1],2) - std::pow(q2i[i],3) / std::pow(q1i[i],2)) + Anorm[i] * (c3 * q3o[i+1] * q2o[i+1] / q1o[i+1] + c4 * std::pow(q2o[i+1],3) * std::pow(q1o[i+1],-2)));
            }

            t += dt;
            swapArrs();
        }

        void step2() {
            for (int i = 0; i < nx+1; i++) {
                int flowLeft = (q2o[i] + q2o[i+1]) < 0;
                
                q1i[i] = q1o[i+flowLeft];
                q2i[i] = q2o[i+flowLeft];
                q3i[i] = q3o[i+flowLeft];
            }

            for (int i = 0; i < nx; i++) {
                q1n[i+1] = q1o[i+1] - dt_dx[i] * (q2i[i+1] - q2i[i] + Anorm[i] * q2o[i+1]);
                q2n[i+1] = q2o[i+1] - dt_dx[i] * (c1 * (q3i[i+1] - q3i[i]) + c2 * (std::pow(q2i[i+1],2)/q1i[i+1] - std::pow(q2i[i],2)/q1i[i]) + Anorm[i] * (c1 * q3o[i+1] + c2 * std::pow(q2o[i+1],2)/q1o[i+1]));
                q3n[i+1] = q3o[i+1] - dt_dx[i] * (c3 * (q3i[i+1]*q2i[i+1]/q1i[i+1] - q3i[i]*q2i[i]/q1i[i]) + c4 * (std::pow(q2i[i+1],3) * std::pow(q1i[i+1],-2) - std::pow(q2i[i],3) * std::pow(q1i[i],-2)) + Anorm[i] * (c3 * q3o[i+1] * q2o[i+1] / q1o[i+1] + c4 * std::pow(q2o[i+1],3) * std::pow(q1o[i+1],-2)));
            }

            t += dt;
            swapArrs();
        }

        void getPressure(float *pArr) {
            for(int i = 0; i < nx; i++) {
                pArr[i] = (q1o[i+1] * q3o[i+1] - std::pow(q2o[i+1],2) / 2) * c1;
            }
        }

        void getVelocity(float *vArr) {
            for (int i = 0; i < nx+1; i++) {
                vArr[i] = q2o[i] / q1o[i];
            }
        }

        float getTime() {
            return t;
        }
};

template <int nx>
class Solver3 {
    private:
        float x[nx + 2];
        float dt_dx[nx];

        float A[nx + 2];
        float Anorm[nx];

        float e[nx+2];

        float q1[2*(nx+2)];
        float q2[2*(nx+2)];

        float q1f[nx+1];
        float q2f[nx+1];
        float ui[nx+1];

        float dt, t, c1;
        float *q1o, *q2o, *q1n, *q2n;

        float ampl = 0.01;
        float freq = 100;

        void swapArrs() {
            bool aOld = q1 == q1o;
            q1o = q1 + (nx+2)*aOld;
            q2o = q2 + (nx+2)*aOld;
            q1n = q1 + (nx+2)*!aOld;
            q2n = q2 + (nx+2)*!aOld;
        }

    public:
        Solver3(float *x, float *A, float *rho, float *rhou, float *e, float dt, float gamma) : dt(dt) {
            q1o = q1; q2o = q2;
            q1n = q1+nx+2; q2n = q2+nx+2;

            c1 = gamma - 1;

            std::copy(x, x+nx, this->x+1);
            std::copy(A, A+nx, this->A+1);
            std::copy(e, e+nx, this->e+1);
            std::copy(rho, rho+nx, this->q1o+1);
            std::copy(rho, rho+nx, this->q1n+1);
            std::copy(rhou, rhou+nx, this->q2o+1);
            std::copy(rhou, rhou+nx, this->q2n+1);
            
            t = 0;

            this->x[0] = 2 * this->x[1] - this->x[2];
            this->x[nx+1] = 2 * this->x[nx] - this->x[nx-1];

            this->A[0] = this->A[1];
            this->A[nx+1] = this->A[nx];

            this->e[0] = this->e[1];
            this->e[nx+1] = this->e[nx];

            for (int i = 0; i < nx; i++) {
                dt_dx[i] = 2 * this->dt / (this->x[i+2] - this->x[i]);
                Anorm[i] = (this->A[i+2] - this->A[i]) / (2 * this->A[i+1]);
            }

            ui[0] = 0;
            q1f[0] = 0;
            q2f[0] = 0;
            ui[nx] = 0;
            q1f[nx] = 0;
            q2f[nx] = 0;
        }

        void setBoundaries(float q1l, float q1r, float q2l, float q2r) {
            q1o[0] = q1l; q1n[0] = q1l;
            q1o[nx+1] = q1r; q1n[nx+1] = q1r;
            q2o[0] = q2l; q2n[0] = q2l;
            q2o[nx+1] = q2r; q2n[nx+1] = q2r;
        }

        void step() {

            for (int i = 0; i < nx+1; i++) {
                ui[i] = (q2o[i] / q1o[i] + q2o[i+1] / q1o[i+1]) / 2;
                int flowLeft = ui[i] < 0;
                
                q1f[i] = q1o[i+flowLeft] * ui[i];
                q2f[i] = std::pow(q2o[i+flowLeft],2) / q1o[i+flowLeft];
            }

            for (int i = 0; i < nx; i++) {
                q1n[i+1] = q1o[i+1] - dt_dx[i] * (q1f[i+1] - q1f[i] + Anorm[i] * q2o[i+1]);
                q2n[i+1] = q2o[i+1] - dt_dx[i] * (q2f[i+1] - q2f[i] + (q1o[i+2] * c1 * e[i+2] - q1o[i] * c1 * e[i]) / 2 + Anorm[i] * (std::pow(q2o[i+1],2)/q1o[i+1] + q1o[i+1] * c1 * e[i+1]));
            }
            
            q1n[0] = q1n[1];
            // q1n[0] = 1.293 + ampl * sin(2 * M_PI * freq * t);
            q1n[nx+1] = q1n[nx];

            // q2n[0] = -q2n[1];
            q2n[0] = -q2n[1];
            q2n[nx+1] = -q2n[nx];

            t += dt;
            swapArrs();
        }

        void getPressure(float *pArr) {
            for (int i = 0; i < nx+2; i++) {
                pArr[i] = q1o[i] * c1 * e[i] * A[i];
            }
        }

        void getVelocity(float *vArr) {
            for (int i = 0; i < nx+1; i++) {
                vArr[i] = q2o[i] / q1o[i];
            }
        }

        float getTime() {
            return t;
        }
};