#pragma once
#include "solver.h"
#include "webster.h"

namespace solver
{

    template <int order>
    void coeff_gen(std::array<float, order> &arr, float r)
    {
        std::array<float, order> coeffs, coeffs_new;
        coeffs.fill(0);
        coeffs_new.fill(0);
        coeffs[0] = 1;

        for (int n = 1; n < order; n += 2)
        {
            for (int i = 0; i < n + 1; i++)
            {
                coeffs_new[n - i] = coeffs[n - i] - r / n * coeffs[i];
            }
            coeffs = coeffs_new;
        }
        arr = coeffs;
    }

    template <int nx, int order>
    class ViscoSolver final : public WebsterIntermediate<ViscoSolver<nx, order>, nx, order + 2>
    {
    private:
        std::array<std::array<float, nx + 1>, order> qvv;
        std::array<std::array<float, nx + 1>, order + 1> qpv;

        std::array<std::array<float, nx + 2>, order> qpp;
        std::array<std::array<float, nx + 2>, order + 1> qvp0;
        std::array<std::array<float, nx + 2>, order + 1> qvp1;
        std::array<float, order + 1> qup;

        friend class Solver<ViscoSolver<nx, order>, nx, 1, 1, order + 2>;

        void step_() override
        {
            std::array<float, nx + 1> &vNext = this->v.getRow(1);

            for (int i = 0; i < nx + 1; i++)
            {
                vNext[i] = -qpv[0][i] * (this->p[i + 1] - this->p[i]);
            }

            for (int r = 1; r < order + 1; r++)
            {
                std::array<float, nx + 1> &v_r = this->v.getRow(-r);
                std::array<float, nx + 2> &p_r = this->p.getRow(-r);
                for (int i = 0; i < nx + 1; i++)
                {
                    vNext[i] -= qvv[r - 1][i] * v_r[i] + qpv[r][i] * (p_r[i + 1] - p_r[i]);
                }
            }

            this->do_bfI();

            this->v.incr();

            std::array<float, nx + 2> &pNext = this->p.getRow(1);

            for (int i = 1; i < nx + 1; i++)
            {
                pNext[i] = -(qvp1[0][i] * this->v[i] - qvp0[0][i] * this->v[i - 1]);
            }

            for (int r = 1; r < order + 1; r++)
            {
                std::array<float, nx + 1> &v_r = this->v.getRow(-r);
                std::array<float, nx + 2> &p_r = this->p.getRow(-r);
                for (int i = 1; i < nx + 1; i++)
                {
                    pNext[i] -= qpp[r - 1][i] * p_r[i] + (qvp1[r][i] * v_r[i] - qvp0[r][i] * v_r[i - 1]);
                }
            }

            this->do_bf();

            this->p.incr();
        }

    public:
        using WebsterIntermediate<ViscoSolver<nx, order>, nx, order + 2>::WebsterIntermediate;
        template <typename T1, typename T2>
        ViscoSolver(
            float k,
            float h,
            float rho,
            float c,
            float eta,
            float nu,
            float gamma,
            std::array<float, nx> &pArr,
            std::array<float, nx + 1> &vArr,
            std::array<float, nx> &SArr,
            T1 &&lbf,
            T2 &&rbf)
            : WebsterIntermediate<ViscoSolver<nx, order>, nx, order + 2>(k, h, pArr, vArr, SArr, lbf, rbf)
        {
            std::cout << "Order: " << order << std::endl;

            // define a and b coefficients for half-derivative approximation
            std::array<float, order + 1> coeff_a, coeff_b;
            coeff_gen<order + 1>(coeff_a, -0.5f);
            coeff_gen<order + 1>(coeff_b, 0.5f);
            std::generate(coeff_b.begin(), coeff_b.end(), [r = 0, coeff_b, k]() mutable
                          { return coeff_b[r++] * std::sqrt(2.0 / k); });

            // declare integer-index LUTs
            std::array<float, nx + 2> a, S, g;

            // declare half-integer-index LUTs
            std::array<float, nx + 1> aI, SI, f, q;

            // generate LUTs for S, a, SI, aI
            std::copy(SArr.begin(), SArr.end(), &(S[1]));
            S[0] = 2 * S[1] - S[2];
            S[nx + 1] = 2 * S[nx] - S[nx - 1];
            std::generate(SI.begin(), SI.end(), [S, i = 0]() mutable
                          { return (S[i] + S[i++ + 1]) / 2; });

            std::generate(a.begin(), a.end(), [S, i = 0]() mutable
                          { return std::sqrt(S[i++] / M_PI); });
            std::generate(aI.begin(), aI.end(), [SI, i = 0]() mutable
                          { return std::sqrt(SI[i++] / M_PI); });

            // generate LUT for g at integer indices
            std::generate(g.begin(), g.end(), [i = 0, gamma, eta, a, nu, c, rho]() mutable
                          { return 2 * (gamma - 1) * std::pow(eta, 0.5) * M_PI * a[i++] / (nu * std::pow(c, 2) * std::pow(rho, 1.5)); });

            // generate LUTs for f, q at half-integer indices
            std::generate(f.begin(), f.end(), [i = 0, rho, eta, aI]() mutable
                          { return 2 * std::pow(rho * eta, 0.5) / aI[i++]; });
            std::generate(q.begin(), q.end(), [i = 0, eta, aI]() mutable
                          { return 3 * eta / std::pow(aI[i++], 2); });

            // generate LUTs for qvv, qpv, qpp, qvp0, qvp1, qup
            float sums[5];

            for (int r = 1; r < order + 1; r++)
            {
                for (int i = 0; i < nx + 1; i++)
                {
                    qvv[r - 1][i] = (coeff_a[r] - coeff_a[r - 1] + k / (2 * rho) * (q[i] * (coeff_a[r] + coeff_a[r - 1]) + f[i] * (coeff_b[r] + coeff_b[r - 1]))) / (1 + k / (2 * rho) * (q[i] + f[i] * coeff_b[0]));
                }
                sums[0] += qvv[r - 1][0];
            }

            for (int r = 0; r < order + 1; r++)
            {
                for (int i = 0; i < nx + 1; i++)
                {
                    qpv[r][i] = ((k * coeff_a[r]) / (rho + k / 2 * (q[i] + f[i] * coeff_b[0]))) / h;
                }
                sums[1] += qpv[r][0];
            }

            for (int r = 1; r < order + 1; r++)
            {
                for (int i = 1; i < nx + 1; i++)
                {
                    qpp[r - 1][i] = (coeff_a[r] - coeff_a[r - 1] + k * rho * std::pow(c, 2) / (2 * S[i]) * g[i] * (coeff_b[r] + coeff_b[r - 1])) / (1 + k * rho * std::pow(c, 2) / (2 * S[i]) * g[i] * coeff_b[0]);
                }
                sums[2] += qpp[r - 1][1];
            }

            for (int r = 0; r < order + 1; r++)
            {
                for (int i = 1; i < nx + 1; i++)
                {
                    float qvp = ((k * rho * std::pow(c, 2) * coeff_a[r]) / (S[i] + k * rho * std::pow(c, 2) / 2 * g[i] * coeff_b[0])) / h;
                    qvp0[r][i] = qvp * SI[i - 1];
                    qvp1[r][i] = qvp * SI[i];
                    qup[r] = qvp * 2;
                }
                sums[3] += qvp0[r][1];
                sums[4] += qvp1[r][1];
            }

            for (int i = 0; i < 5; i++)
            {
                std::cout << sums[i] << ", ";
            }
            std::cout << std::endl;
        }

        
    };

    template <int nx, int order>
    struct BFuncOpen<ViscoSolver<nx, order>> : public BFuncWIOpen<ViscoSolver<nx, order>> {};

    template <int nx, int order>
    struct BFuncClosed<ViscoSolver<nx, order>> : public BFuncWIClosed<ViscoSolver<nx, order>> {};

    template <int nx, int order>
    struct BFuncDriver<ViscoSolver<nx, order>> : public BFuncWIDriver<ViscoSolver<nx, order>> {};
}