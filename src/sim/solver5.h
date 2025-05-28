#include <math.h> 
#include <algorithm>
#include <iostream>
#include <functional>
#include <vector>
#include <chrono>

enum BoundaryType {
    Zero, Closed, Open, In, Out, Driver
};

template <int nx, int m>
class Solver5 {
private:

    float h, k;

    float qvv[m][nx+1], qpv[m+1][nx+1], qpp[m][nx+2], qvp0[m+1][nx+2], qvp1[m+1][nx+2], qup[m+1];

    // pressure array
    float p[(nx+2) * (m+2)];
    // velocity array
    float v[(nx+1) * (m+2)];

    float * p_curr;
    float * p_next;

    float * v_curr;
    float * v_next;

    long step_count = 0;

    BoundaryType lb, rb;

    float driver_freq = 440;
    float driver_ampl = 1e1;

    void (Solver5<nx, m>::*lbfp)(float *pi, float *po),
        (Solver5<nx, m>::*rbfp)(float *pi, float *po),
        (Solver5<nx, m>::*lbfv)(float *vi, float *vo),
        (Solver5<nx, m>::*rbfv)(float *vi, float *vo);

    void coeff_gen(float * arr, float r, int order) {
        float coeffs[order], coeffs_new[order];
        std::fill(coeffs, coeffs+order, 0);
        std::fill(coeffs_new, coeffs_new+order, 0);
        coeffs[0] = 1;
        for(int n = 1; n < order; n += 2) {
            for(int i = 0; i < n+1; i++) {
                coeffs_new[n-i] = coeffs[n-i] - r/n*coeffs[i];
            }
            std::copy(coeffs_new, coeffs_new+order, coeffs);
        }
        std::copy(coeffs, coeffs + order, arr);
    }

    // move circular queue head pointers
    void step_cq_p() {
        p_curr = p_next;
        p_next = (p_next - p + nx+2) % ((nx+2) * (m+2)) + p;
    }

    void step_cq_v() {
        v_curr = v_next;
        v_next = (v_next - v + nx+1) % ((nx+1) * (m+2)) + v;
    }

    // get pressure pointer from n steps prior
    float * get_prev_p_ptr(int n) {
        return (p_curr - p - n*(nx+2) + (nx+2) * (m+2)) % ((nx+2) * (m+2)) + p;
    }

    float * get_prev_v_ptr(int n) {
        return (v_curr - v - n*(nx+1) + (nx+1) * (m+2)) % ((nx+1) * (m+2)) + v;
    }

    void bp_equate(float *pi, float *po) {*po = *pi;}

    void bp_sin_driver(float *pi, float *po) {
        /*
        po = -qup[0] * v;
        for (int r = 1; r < M+1; r++) {
            po -= qpp[r-1][0] - qup[r];
        }
        */
        *po = *pi;
    }

    void bp_zero(float *pi, float *po) {*po = 0;}

    void bv_zero(float *vi, float *vo) {*vo = 0;}

    void bv_pass(float *vi, float *vo) {}

    void bv_sin_driver(float *vi, float *vo) {
        *vo = driver_ampl * std::sin(2*M_PI*driver_freq*step_count*k);
    }

    // enforce boundary functions
    void enforce_bfp() {
        (this->*lbfp)(p_curr+1, p_curr);
        (this->*rbfp)(p_curr+nx, p_curr+nx+1);
    }

    void enforce_bfv() {
        (this->*lbfv)(v_curr+1, v_curr);
        (this->*rbfv)(v_curr+nx-1, v_curr+nx);
    }

public:

    Solver5(float h, float k, float c, float rho, float eta, float nu, float gamma, float * p, float * v, float * a_in, BoundaryType lb, BoundaryType rb) : step_count(0), h(h), k(k), lb(lb), rb(rb) {
        
        std::cout << "CFL: " << c*k/h << std::endl;
        std::cout << "Order: " << m << std::endl;

        // define a and b coefficients for half-derivative approximation
        float coeff_a[m+1], coeff_b[m+1];
        coeff_gen(coeff_a, -0.5, m+1);
        coeff_gen(coeff_b, 0.5, m+1);
        std::generate(coeff_b, coeff_b+m+1, [r=0, coeff_b, k]() mutable {return coeff_b[r++] * std::sqrt(2.0/k);});

        float a[nx+2], S[nx+2], g[nx+2], ai[nx+1], Si[nx+1], f[nx+1], q[nx+1];

        // generate lookup tables for a, g at integer indices
        std::copy(a_in, a_in+nx, a+1);
        a[0] = 2 * a[1] - a[2];
        a[nx+1] = 2 * a[nx] - a[nx-1];

        std::generate(S, S+nx+1+2, [i=0, a]() mutable {return M_PI * std::pow(a[i++], 2);});
        
        std::generate(g, g+nx+2, [i=0, gamma, eta, a, nu, c, rho]() mutable {return 2*(gamma-1)*std::pow(eta,0.5)*M_PI*a[i++]/(nu*std::pow(c,2)*std::pow(rho,1.5));});

        // generate lookup tables for a, f, q at half-integer indices
        std::generate(ai, ai+nx+1, [i=0,a]() mutable {return 0.5 * (a[i] + a[i++ +1]);});
        std::generate(Si, Si+nx+1, [i=0, S]() mutable {return 0.5 * (S[i] + S[i++ +1]);});
        std::generate(f, f+nx+1, [i=0, rho, eta, ai]() mutable {return 2*std::pow(rho*eta,0.5)/ai[i++];});
        std::generate(q, q+nx+1, [i=0, eta, ai]() mutable {return 3*eta/std::pow(ai[i++],2);});
        
        float sums[5];

        for (int r = 1; r < m+1; r++) {
            for (int i = 0; i < nx+1; i++) {
                qvv[r-1][i] = (coeff_a[r] - coeff_a[r-1] + k/(2*rho)*(q[i]*(coeff_a[r]+coeff_a[r-1])+f[i]*(coeff_b[r]+coeff_b[r-1]))) / (1 + k/(2*rho)*(q[i]+f[i]*coeff_b[0]));
            }
            sums[0] += qvv[r-1][0];
        }

        for (int r = 0; r < m+1; r++) {
            for (int i = 0; i < nx+1; i++) {
                qpv[r][i] = ((k * coeff_a[r]) / (rho + k/2*(q[i]+f[i]*coeff_b[0]))) / h;
            }
            sums[1] += qpv[r][0];
        }

        for (int r = 1; r < m+1; r++) {
            for (int i = 1; i < nx+1; i++) {
                qpp[r-1][i] = (coeff_a[r] - coeff_a[r-1] + k*rho*std::pow(c,2)/(2*S[i])*g[i]*(coeff_b[r]+coeff_b[r-1])) / (1+k*rho*std::pow(c,2)/(2*S[i])*g[i]*coeff_b[0]);
            }
            sums[2] += qpp[r-1][1];
        }

        for (int r = 0; r < m+1; r++) {
            for (int i = 1; i < nx+1; i++) {
                float qvp = ((k*rho*std::pow(c,2)*coeff_a[r]) / (S[i]+k*rho*std::pow(c,2)/2*g[i]*coeff_b[0])) / h;
                qvp0[r][i] = qvp * Si[i-1];
                qvp1[r][i] = qvp * Si[i];
                qup[r] = qvp * 2;
            }
            sums[3] += qvp0[r][1];
            sums[4] += qvp1[r][1];
        }

        for (int i = 0; i < 5; i++) {
            std::cout << sums[i] << ",";
        }
        std::cout << std::endl;

        p_curr = this->p;
        v_curr = this->v;
        p_next = (p_curr + nx+2 - this->p) % ((nx+2) * (m+2)) + this->p;
        v_next = (v_curr + nx+1 - this->v) % ((nx+1) * (m+2)) + this->v;

        
        for (int i = 0; i < m + 1; i++) {
            std::copy(p, p+nx, this->p+i*(nx+2)+1);
            std::copy(v, v+nx-1, this->v+i*(nx+1)+1);

            this->p[i*(nx+2)] = 2 * this->p[1+i*(nx+2)] - this->p[2+i*(nx+2)];
            this->p[nx+1+i*(nx+2)] = 2 * this->p[nx+i*(nx+2)] - this->p[nx-1+i*(nx+2)];

            this->v[i*(nx+1)] = 2 * this->v[1+i*(nx+1)] - this->v[2+i*(nx+1)];
            this->v[nx+i*(nx+1)] = 2 * this->v[nx-1+i*(nx+1)] - this->v[nx-2+i*(nx+1)];
        }

        switch(this->lb) {
            case BoundaryType::Closed:
                lbfp = bp_equate;
                lbfv = bv_zero;
                break;
            case BoundaryType::Open:
                lbfp = bp_equate;
                lbfv = bv_pass;
                break;
            case BoundaryType::Driver:
                lbfp = bp_sin_driver;
                lbfv = bv_sin_driver;
                break;
            case BoundaryType::Zero:
                lbfp = bp_zero;
                lbfv = bv_zero;
                break;
            default:
                #pragma message("Boundary condition not implemented!")
                lbfp = bp_equate;
                lbfv = bv_pass;
        }

        switch(this->rb) {
            case BoundaryType::Closed:
                rbfp = bp_equate;
                rbfv = bv_zero;
                break;
            case BoundaryType::Open:
                rbfp = bp_equate;
                rbfv = bv_pass;
                break;
            case BoundaryType::Driver:
                rbfp = bp_sin_driver;
                rbfv = bv_sin_driver;
                break;
            case BoundaryType::Zero:
                rbfp = bp_zero;
                rbfv = bv_zero;
                break;
            default:
            #pragma message("Boundary condition not implemented!")
                rbfp = bp_equate;
                rbfv = bv_pass;
        }
    }

    void step() {
        
        auto start = std::chrono::high_resolution_clock::now();
        enforce_bfv();
        for (int i = 0; i < nx+1; i++) {
            v_next[i] = -qpv[0][i]*(p_curr[i+1]-p_curr[i]);
        }

        for(int r = 1; r < m + 1; r++) {
            float *prev_v_ptr = get_prev_v_ptr(r);
            float *prev_p_ptr = get_prev_p_ptr(r);

            for (int i = 0; i < nx+1; i++) {
                v_next[i] -= qvv[r-1][i]*prev_v_ptr[i] + qpv[r][i]*(prev_p_ptr[i+1] - prev_p_ptr[i]);
            }
        }
        step_cq_v();
        auto end = std::chrono::high_resolution_clock::now();
        double timeTaken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * 1e-9;
        
        enforce_bfp();
        for (int i = 1; i < nx+1; i++) {
            p_next[i] = -(qvp1[0][i]*v_curr[i] - qvp0[0][i]*v_curr[i-1]);
        }

        for(int r = 1; r < m + 1; r++) {
            float *prev_v_ptr = get_prev_v_ptr(r);
            float *prev_p_ptr = get_prev_p_ptr(r);

            for (int i = 1; i < nx+1; i++) {
                p_next[i] -= qpp[r-1][i]*prev_p_ptr[i] + (qvp1[r][i]*prev_v_ptr[i] - qvp0[r][i]*prev_v_ptr[i-1]);
            }
        }
        step_cq_p();
        step_count++;
    
    }

    void getPressure(float * pArr) {
        std::copy(p_curr, p_curr+nx+2, pArr);
    }

    void getVelocity(float * vArr) {
        std::copy(v_curr, v_curr+nx+1, vArr);
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
        return v_curr[i];
    }

    float getTime() {
        return step_count * k;
    }
};