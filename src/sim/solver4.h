#include <math.h> 
#include <algorithm>
#include <iostream>
#include <functional>
#include <vector>

enum BoundaryType {
    Zero, Closed, Open, In, Out, Driver
};

template <int nx>
class Solver4 {
private:

    float h, k, c_squared, rho;

    // pressure array
    float p[(nx+2) * (2)];
    float v[nx+1];
    float S[nx+2];
    float Si[nx+1];

    float * p_curr;
    float * p_next;

    long step_count;

    float driver_freq = 220;
    float driver_ampl = 10;

    BoundaryType lb, rb;

    void (Solver4<nx>::*lbf)(float * pi, float * po, float * v), (Solver4<nx>::*rbf)(float * pi, float * po, float * v);

    // move circular queue head pointer for pressure array
    void step_cq() {
        p_curr = p_next;
        p_next = (p_next + nx+2 - p) % ((nx+2) * (2)) + p;
    }

    // get pressure pointer from n steps prior
    float * get_prev_p_ptr(int n) {
        return (p_curr - n*(nx+2) - p) % ((nx+2) * (2)) + p;
    }

    // closed boundary condition
    void b_closed(float * pi, float * po, float * v) {
        *po = *pi;
        *v = 0;
    }

    // open boundary condition
    void b_open(float * pi, float * po, float * v) {
        *po = *pi;
    }

    // enforce boundary functions
    void enforce_bfs() {
        (this->*lbf)(p_next+1, p_next, v);
        (this->*rbf)(p_next+nx, p_next+nx+1, v+nx);
    }

    void b_driver(float * pi, float * po, float * v) {
        *po = *pi;
        *v = driver_ampl * std::sin(2*M_PI*driver_freq*step_count*k);
    }

public:

    Solver4(float h, float k, float c, float rho, float * p, float * v, float * S, BoundaryType lb, BoundaryType rb) : h(h), k(k), c_squared(std::pow(c,2)), rho(rho), lb(lb), rb(rb), step_count(0) {
        p_curr = this->p;
        p_next = (p_curr + nx+2 - this->p) % ((nx+2) * (2)) + this->p;
        
        std::copy(p, p+nx, p_curr+1);
        std::copy(v, v+nx-1, this->v+1);
        std::copy(S, S+nx, this->S+1);
        
        p_curr[0] = 2 * p_curr[1] - p_curr[2];
        p_curr[nx+1] = 2 * p_curr[nx] - p_curr[nx-1];

        this->v[0] = 2 * this->v[1] - this->v[2];
        this->v[nx] = 2 * this->v[nx-1] - this->v[nx-2];

        this->S[0] = 2 * this->S[1] - this->S[2];
        this->S[nx+1] = 2 * this->S[nx] - this->S[nx-1];

        for (int i = 0; i < nx+1; i++) {
            Si[i] = 0.5 * (S[i] + S[i+1]);
        }

        switch(this->lb) {
            case BoundaryType::Closed:
                lbf = b_closed;
                break;
            case BoundaryType::Open:
                lbf = b_open;
                break;
            case BoundaryType::Driver:
                lbf = b_driver;
                break;
            default:
                #pragma message(": warning<b_cond_not_implemented> Boundary condition not implemented!")
                lbf = b_open;
        }

        switch(this->rb) {
            case BoundaryType::Closed:
                rbf = b_closed;
                break;
            case BoundaryType::Open:
                rbf = b_open;
                break;
            case BoundaryType::Driver:
                rbf = b_driver;
                break;
            default:
                #pragma message(": warning<b_cond_not_implemented> Boundary condition not implemented!")
                rbf = b_open;
        }
    }

    void step() {
        for (int i = 0; i < nx+1; i++) {
            v[i] -= k / (h * rho) * (p_curr[i+1] - p_curr[i]);
        }

        for (int i = 0; i < nx; i++) {
            p_next[i+1] = p_curr[i+1] - rho * c_squared * k / (h * S[i+1]) * (Si[i+1] * v[i+1] - Si[i] * v[i]);
        }

        enforce_bfs();
        step_count++;
        step_cq();
    }

    void getPressure(float * pArr) {
        std::copy(p_curr, p_curr+nx+2, pArr);
    }

    void getVelocity(float * vArr) {
        std::copy(v, v+nx+1, vArr);
    }

    float getPressure(int i) {
        if (i < 0 || nx+2 < i) {
            return NAN;
        }
        return p_curr[i];
    }

    float getVelocity(int i) {
        if (i < 0 || nx+1 < i) {
            return NAN;
        }
        return v[i];
    }

    float getTime() {
        return step_count * k;
    }
};