#include <math.h>

template <int nx>
class Solver {
    private:
        float x[nx];
        float xi[nx+1];
        float rho[nx];
        float rhou[nx];
        float e[nx];

        float gamma;
        float time;

        float ui[nx+1];
        float fluxrho[nx+1];
        float fluxrhou[nx+1];
        float p[nx];

    public:
        Solver(float *x, float *rho, float *rhou, float *e, float gamma) : gamma(gamma) {
            std::copy(x, x+nx, this->x);
            std::copy(rho, rho+nx, this->rho);
            std::copy(rhou, rhou+nx, this->rhou);
            std::copy(e, e+nx, this->e);

            time = 0;
            for (int i = 1; i < nx; i++) {
                xi[i] = 0.5 * (x[i]+x[i-1]);
            }
            xi[0] = 2*xi[1] - xi[2];
            xi[nx] = 2*xi[nx-1] - xi[nx-2];
        }

        void step(float dt) {
            
            // compute cell interface velocities as mean of adjacent cell velocities
            for (int i = 1; i < nx; i++) {
                ui[i] = 0.5 * (rhou[i] / rho[i] + rhou[i-1] / rho[i-1]);
            }

            // compute density flux
            for (int i = 1; i < nx; i++) {
                fluxrho[i] = rho[i-(ui[i]>0)]*ui[i];
            }

            // update density
            for (int i = 0; i < nx; i++) {
                rho[i] -= dt/(xi[i+1]-xi[i]) * (fluxrho[i+1]-fluxrho[i]);
            }

            // compute momentum density flux
            for (int i = 1; i < nx; i++) {
                fluxrhou[i] = std::pow(rhou[i-(ui[i]>0)], 2) / rho[i-(ui[i]>0)];
            }

            // update momentum density
            for (int i = 0; i < nx; i++) {
                rhou[i] -= dt/(xi[i+1]-xi[i]) * (fluxrhou[i+1]-fluxrhou[i]);
            }

            for (int i = 0; i < nx; i++) {
                p[i] = rho[i] * e[i] * (gamma - 1);
            }

            for (int i = 1; i < nx-1; i++) {
                rhou[i] -= dt * (p[i+1]-p[i-1])/(x[i+1]-x[i-1]);
            }

            rhou[0] -= 0.5*dt*(p[1]-p[0])/(x[1]-x[0]);
            rhou[nx-1] -= 0.5*dt*(p[nx-1]-p[nx-2])/(x[nx-1]-x[nx-2]);

            time += dt;
        }

        void getPressure(float *pArr) {
            std::copy(p, p+nx, pArr);
        }

        float getTime() {
            return time;
        }
};