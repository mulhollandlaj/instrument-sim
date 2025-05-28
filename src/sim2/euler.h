#pragma once
#include "solver.h"
namespace solver
{
    template <int nx>
    class EulerSolver : public Solver<EulerSolver<nx>, nx, 2, 0, 2>
    {
    protected:
        std::array<float, nx + 2> S;
        std::array<float, nx> SNorm;
        VarArr<float, nx + 2, 2> &rho;
        VarArr<float, nx + 2, 2> &rho_u;

        std::array<float, nx + 2> e;

        std::array<float, nx + 1> rhoFlux;
        std::array<float, nx + 1> rho_uFlux;

        const float dt_dx, c1, p0, rho0;

        friend Solver<EulerSolver<nx>, nx, 2, 0, 1>;
        friend BFuncOpen<EulerSolver<nx>>;
        friend BFuncClosed<EulerSolver<nx>>;
        friend BFuncDriver<EulerSolver<nx>>;

        void step_() override {
            for (int i = 0; i < nx+1; i++) {
                float uI = (rho_u[i] / rho[i] + rho_u[i+1] / rho[i+1]) / 2;
                int flowLeft = uI < 0;
                
                rhoFlux[i] = rho[i+flowLeft] * uI;
                rho_uFlux[i] = std::pow(rho_u[i+flowLeft],2) / rho[i+flowLeft];
            }

            for (int i = 0; i < nx; i++) {
                rho.getRow(1)[i+1] = rho[i+1] - dt_dx * (rhoFlux[i+1] - rhoFlux[i] + SNorm[i] * rho_u[i+1]);
                rho_u.getRow(1)[i+1] = rho_u[i+1] - dt_dx * (rho_uFlux[i+1] - rho_uFlux[i] + (rho[i+2] * c1 * e[i+2] - rho[i] * c1 * e[i]) / 2 + SNorm[i] * (std::pow(rho_u[i+1],2)/rho[i+1] + rho[i+1] * c1 * e[i+1]));
            }

            this->do_bf();

            rho.incr();
            rho_u.incr();
        }

    public:
        using Solver<EulerSolver<nx>, nx, 2, 0, 2>::Solver;
        template <typename T1, typename T2>
        EulerSolver(
            float k_,
            float h_,
            float gamma,
            float p0_,
            float rho0_,
            std::array<float, nx> &rho_,
            std::array<float, nx> &rho_u_,
            std::array<float, nx> &S_,
            std::array<float, nx> &e_, 
            T1 &&lbf_,
            T2 &&rbf_)
            : Solver<EulerSolver<nx>, nx, 2, 0, 2>::Solver(
                  k_,
                  h_,
                  std::array<std::array<float, nx>, 2>{rho_, rho_u_},
                  std::array<std::array<float, nx + 1>, 0>{},
                  std::make_shared<std::decay_t<T1>>(std::forward<T1>(lbf_)),
                  std::make_shared<std::decay_t<T2>>(std::forward<T2>(rbf_))),
              rho(this->q[0]),
              rho_u(this->q[1]),
              dt_dx(this->k / this->h),
              c1(gamma - 1),
              p0(p0_),
              rho0(rho0)
        {
            std::copy(S_.begin(), S_.end(), S.begin() + 1);
            S[0] = 2 * S[1] - S[2];
            S[nx+1] = 2 * S[nx] - S[nx-1];
            std::generate(SNorm.begin(), SNorm.end(), [this, i = 0]() mutable
                          { return (this->S[i + 2] - this->S[i]) / (2 * this->S[i++ + 1]); });

            std::copy(e_.begin(), e_.end(), e.begin() + 1);
            e[0] = 2 * e[1] - e[2];
            e[nx+1] = 2 * e[nx] - e[nx-1];

            rho[0] = rho[1];
            rho[nx+1] = rho[nx];

            rho_u[0] = rho_u[1];
            rho_u[nx+1] = rho_u[nx];

            std::generate(rho_.begin(), rho_.end(), [rho_, this, i = 0]() mutable {return (rho_[i] / this->S[i++]);});
            std::generate(rho_u_.begin(), rho_u_.end(), [rho_u_, this, i = 0]() mutable {return (rho_u_[i] / this->S[i++]);});
        }

        void getPressure(std::array<float, nx> &destArr) const override
        {
            std::generate(destArr.begin(), destArr.end(), [this, i = 1]() mutable
                          { return (this->rho[i] * c1 * e[i] * S[i++]) - p0; });
        }

        void getVelocity(std::array<float, nx> &destArr) const override
        {
            std::generate(destArr.begin(), destArr.end(), [this, i = 1]() mutable
                          { return this->rho_u[i] / this->rho[i] * S[i++]; });
        }

        float getPressure(float x) const override
        {
            return this->arrayLerp(x, rho.getRow(0)) * c1 * this->arrayLerp(x, e) * this->arrayLerp(x, S) - p0;
        }

        float getVelocity(float x) const override
        {
            return this->arrayLerp(x, rho_u.getRow(0)) / this->arrayLerp(x, rho.getRow(0)) * this->arrayLerp(x, S);
        }
    };

    template <int nx>
    struct BFuncOpen<EulerSolver<nx>> : public BFunc<EulerSolver<nx>> {
        void bf(EulerSolver<nx> &s, int inIx, int outIx)
        {
            s.rho.getRow(1)[outIx] = s.rho0;
            s.rho_u.getRow(1)[outIx] = s.rho_u.getRow(1)[inIx];
        }

        void bfI(EulerSolver<nx> &s, int inIx, int outIx) {}
    };
    
    template <int nx>
    struct BFuncClosed<EulerSolver<nx>> : public BFunc<EulerSolver<nx>> {
        void bf(EulerSolver<nx> &s, int inIx, int outIx)
        {
            s.rho.getRow(1)[outIx] = s.rho.getRow(1)[inIx];
            s.rho_u.getRow(1)[outIx] = -s.rho_u.getRow(1)[inIx];
        }
        
        void bfI(EulerSolver<nx> &s, int inIx, int outIx) {}
    };
    
        template <int nx>
        struct BFuncDriver<EulerSolver<nx>> : public Phasor, public BFunc<EulerSolver<nx>> {
            void bf(EulerSolver<nx> &s, int inIx, int outIx)
            {
                s.rho.getRow(1)[outIx] = s.rho.getRow(1)[inIx];
                s.rho_u.getRow(1)[outIx] = s.rho.getRow(1)[outIx] * getValue(s.getTime());
            }
    
            void bfI(EulerSolver<nx> &s, int inIx, int outIx) {}

            BFuncDriver<EulerSolver<nx>>(float ampl_, float freq_, float phase_) : Phasor(ampl_, freq_, phase_) {}
        };
}