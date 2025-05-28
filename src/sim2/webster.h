#pragma once
#include "solver.h"

namespace solver
{
    template <typename Derived>
    struct BFuncWIOpen;

    template <typename Derived>
    struct BFuncWIClosed;

    template <typename Derived>
    struct BFuncWIDriver;

    template <typename Derived>
    struct is_WI_derived;

    template <typename D, int nx,  int r>
    class WebsterIntermediate : public Solver<D, nx, 1, 1, r>
    {
    protected:
        std::array<float, nx + 2> S;
        std::array<float, nx + 1> SI;
        VarArr<float, nx + 2, r> &p;
        solver::VarArr<float, nx + 1, r> &v;

        friend Solver<D, nx, 1, 1, r>;
        friend BFuncWIOpen<D>;
        friend BFuncWIClosed<D>;
        friend BFuncWIDriver<D>;

    public:
        using Solver<D, nx, 1, 1, r>::Solver;
        template <typename T1, typename T2>
        WebsterIntermediate(
            float k_,
            float h_,
            std::array<float, nx> &pArr,
            std::array<float, nx + 1> &vArr,
            std::array<float, nx> &SArr,
            T1 &&lbf_,
            T2 &&rbf_)
            : Solver<D, nx, 1, 1, r>(
                  k_,
                  h_,
                  std::array<std::array<float, nx>, 1>{pArr},
                  std::array<std::array<float, nx + 1>, 1>{vArr},
                  std::make_shared<std::decay_t<T1>>(std::forward<T1>(lbf_)),
                  std::make_shared<std::decay_t<T2>>(std::forward<T2>(rbf_))),
              p(this->q[0]),
              v(this->qI[0])
        {
            static_assert(std::is_base_of<WebsterIntermediate, D>::value,
                          "CRTP misused: Derived must inherit from WebsterIntermediate<Derived, ...>");
            std::copy(SArr.begin(), SArr.end(), &(S[1]));
            S[0] = S[1];
            S[nx + 1] = S[nx];
            std::generate(&(SI[0]), &(SI[nx + 1]), [this, i = 0]() mutable
                          { return (S[i] + S[i++ + 1]) / 2; });
        }

        virtual void getPressure(std::array<float, nx> &destArr) const override
        {
            std::copy(&(p[1]), &(p[nx + 1]), destArr.begin());
        }

        virtual void getVelocity(std::array<float, nx> &destArr) const override
        {
            std::generate(destArr.begin(), destArr.end(), [this, i = 0]() mutable
                          { return (this->v[i] + this->v[i++ + 1]) / 2; });
        }

        virtual float getPressure(float x) const override
        {
            return this->arrayLerp(x, p.getRow(0));
        }

        virtual float getVelocity(float x) const override
        {
            return this->arrayLerpI(x, v.getRow(0));
        }
    };

    template <typename Derived>
    struct is_WI_derived
    {
    private:
        template <int nx, int r>
        static std::true_type test(const WebsterIntermediate<Derived, nx, r> *);

        static std::false_type test(...);

    public:
        static constexpr bool value = decltype(test(std::declval<Derived *>()))::value;
    };

    template <typename Derived>
    struct BFuncWI : public virtual BFunc<Derived>
    {
        BFuncWI()
        {
            static_assert(is_WI_derived<Derived>::value, "Derived must derive from WebsterIntermediate");
        }
    };

    template <typename Derived>
    struct BFuncWIOpen : public BFuncWI<Derived>
    {
        void bf(Derived &solver, int inIx, int outIx) override
        {
            solver.p.getRow(1)[outIx] = 0;
        }

        void bfI(Derived &solver, int inIx, int outIx) override
        {
            // solver.v.getRow(1)[outIx] = solver.v.getRow(1)[inIx];
        }
    };

    template <typename Derived>
    struct BFuncWIClosed : public BFuncWI<Derived>
    {
        void bf(Derived &solver, int inIx, int outIx) override
        {
            solver.p.getRow(1)[outIx] = solver.p.getRow(1)[inIx];
        }

        void bfI(Derived &solver, int inIx, int outIx) override
        {
            solver.v.getRow(1)[outIx] = 0;
        }
    };

    template <typename Derived>
    struct BFuncWIDriver : public BFuncWI<Derived>, public Phasor
    {
        void bf(Derived &solver, int inIx, int outIx) override
        {
            solver.p.getRow(1)[outIx] = solver.p.getRow(1)[inIx];
        }

        void bfI(Derived &solver, int inIx, int outIx) override
        {
            solver.v.getRow(1)[outIx] = getValue(solver.getTime());
        }
    };

    template <int nx>
    class WebsterSolver final : public WebsterIntermediate<WebsterSolver<nx>, nx, 2>
    {
    private:
        float rho, c_squared;

        void step_() override
        {
            for (int i = 0; i < nx + 1; i++)
            {
                this->v.getRow(1)[i] = this->v[i] - this->k / (this->h * rho) * (this->p[i + 1] - this->p[i]);
            }
            this->do_bfI();
            this->v.incr();

            for (int i = 0; i < nx; i++)
            {
                this->p.getRow(1)[i + 1] = this->p[i + 1] - rho * c_squared * this->k / (this->h * this->S[i + 1]) * (this->SI[i + 1] * this->v[i + 1] - this->SI[i] * this->v[i]);
            }
            this->do_bf();
            this->p.incr();
        }

    public:
        using WebsterIntermediate<WebsterSolver<nx>, nx, 2>::WebsterIntermediate;

        template <typename T1, typename T2>
        WebsterSolver(
            float k_,
            float h_,
            float rho_,
            float c_,
            std::array<float, nx> &pArr,
            std::array<float, nx + 1> &vArr,
            std::array<float, nx> &SArr,
            T1 &&lbf_,
            T2 &&rbf_)
            : WebsterIntermediate<WebsterSolver<nx>, nx, 2>(
                  k_,
                  h_,
                  pArr,
                  vArr,
                  SArr,
                  std::forward<T1 &&>(lbf_),
                  std::forward<T2 &&>(rbf_)),
              rho(rho_),
              c_squared(std::pow(c_, 2)) {}
    };

    template <int nx>
    struct BFuncOpen<WebsterSolver<nx>> : public BFuncWIOpen<WebsterSolver<nx>>
    {
    };

    template <int nx>
    struct BFuncClosed<WebsterSolver<nx>> : public BFuncWIClosed<WebsterSolver<nx>>
    {
    };

    template <int nx>
    struct BFuncDriver<WebsterSolver<nx>> : public BFuncWIDriver<WebsterSolver<nx>>
    {
    };
}