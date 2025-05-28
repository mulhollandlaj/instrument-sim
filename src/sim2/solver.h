#pragma once
#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <sstream>
#include <array>
#include <utility>
#include <unordered_map>
#include <memory>
#include <iostream>

namespace solver
{
    template <typename T, int rowSize, int rows>
    class VarArr
    {
    private:
        std::array<std::array<T, rowSize>, rows> arr;
        int h;

    public:
        VarArr()
        {
            h = 0;
        }

        void incr()
        {
            h = ((h + 1 + rows) % rows);
        }

        inline std::array<T, rowSize> &getRow(int rowOffset)
        {
#ifdef DEBUG
            if (rowOffset < 2 - rows)
            {
                std::cout << "VarArr accessed illegal row" << std::endl;
            }
#endif
            return arr[(h + rows + rowOffset) % rows];
        }

        T &get(int rowOffset, int ix)
        {
            return arr[(h + rows + rowOffset) % rows][ix];
        }

        const T &get(int rowOffset, int ix) const
        {
            return arr[(h + rows + rowOffset) % rows][ix];
        }

        T &operator[](int ix)
        {
            return arr[h][ix];
        }

        const T &operator[](int ix) const
        {
            return arr[h][ix];
        }
    };

    class Phasor
    {
    private:
        float freq;
        float ampl;
        float phase;

    public:
        Phasor(float freq_, float ampl_, float phase_) : freq(freq_), ampl(ampl_), phase(phase_) {}

        float getFreq() const { return freq; }
        float getAmpl() const { return ampl; }
        float getPhase() const { return phase; }

        void setFreq(float freq_) { freq = freq_; }
        void setAmpl(float ampl_) { ampl = ampl_; }
        void setPhase(float phase_) { phase = fmod(fmod(phase_, (2 * M_PI)) + 2 * M_PI, 2 * M_PI); }

        float getValue(float t) const
        {
            return ampl * std::sin(2 * M_PI * freq * t + phase);
        }
    };

    template <class D>
    struct BFunc;

    template <class D>
    struct BFuncOpen;

    template <class D>
    struct BFuncClosed;

    template <class D>
    struct BFuncDriver;

    template <class D, int nx, int nv, int ni, int r>
    class Solver
    {
    friend BFuncOpen<D>;
    friend BFuncClosed<D>;
    friend BFuncDriver<D>;

    private:
        D &asD()
        {
            return static_cast<D &>(*this);
        }

    protected:
        const float k, h;
        long ns = 0;

        // std::array<float, nx + 2> x;
        // std::array<float, nx + 1> xI;

        std::array<VarArr<float, nx + 2, r>, nv> q;
        std::array<VarArr<float, nx + 1, r>, ni> qI;

        std::shared_ptr<BFunc<D>> lbf, rbf;

        
        void do_bf()
        {
            lbf->bf(asD(), 1, 0);
            rbf->bf(asD(), nx, nx + 1);
        }

        void do_bfI()
        {
            lbf->bfI(asD(), 1, 0);
            rbf->bfI(asD(), nx - 1, nx);
        }

        float arrayLerp(float x, const std::array<float, nx + 2> &arr) const
        {
            float floatIx = 1 + x / h;
            int ix = (int)floatIx;
            if (ix > nx)
            {
                return arr[nx+1];
            }

            if (ix < 1)
            {
                return arr[0];
            }

            float offset = floatIx - ix;

            return offset * arr[ix + 1] + (1 - offset) * arr[ix];
        }

        float arrayLerpI(float x, const std::array<float, nx + 1> &arr) const
        {
            float floatIx = 0.5 + x / h;
            int ix = (int)floatIx;
            if (ix > nx - 1 || ix < 0)
            {
                return NAN;
            }
            float offset = floatIx - ix;

            return offset * arr[ix + 1] + (1 - offset) * arr[ix];
        }

        virtual void step_() = 0;

    public:
        ~Solver() = default;
        explicit Solver(
            float k_,
            float h_,
            std::array<std::array<float, nx>, nv> qArrs,
            std::array<std::array<float, nx + 1>, ni> qArrsI,
            std::shared_ptr<BFunc<D>> lbf_, std::shared_ptr<BFunc<D>> rbf_)
            : k(k_), h(h_), lbf(std::move(lbf_)), rbf(std::move(rbf_))
        {
            static_assert(std::is_base_of<Solver, D>::value,
                          "CRTP misused: Derived must inherit from Solver<Derived, ...>");
/*
            std::generate(&(x[0]), &(x[nx + 2]), [this, i = -1]() mutable
                          { return this->h * i++; });
            std::generate(&(xI[0]), &(xI[nx + 1]), [this, i = 0]() mutable
                          { return this->h * (i++ - 0.5); });
*/
            for (int i = 0; i < nv; i++)
            {
                std::copy(&(qArrs[i][0]), &(qArrs[i][nx]), &(q[i][1]));
            }

            for (int i = 0; i < ni; i++)
            {
                std::copy(&(qArrsI[i][0]), &(qArrsI[i][nx + 1]), &(qI[i][0]));
            }
        }

        float getTime() const
        {
            return ns * k;
        }

        VarArr<float, nx + 2, r> &getVar(int varIndex)
        {
            return q[varIndex];
        }

        VarArr<float, nx + 1, r> &getVarI(int varIndex)
        {
            return qI[varIndex];
        }

        BFunc<D> &getBFunc(bool right) const
        {
            return right ? rbf : lbf;
        }

        void step()
        {
            step_();
            ns++;
        }

        virtual void getPressure(std::array<float, nx> &destArr) const = 0;
        virtual void getVelocity(std::array<float, nx> &destArr) const = 0;
        virtual float getPressure(float x) const = 0;
        virtual float getVelocity(float x) const = 0;
    };

    template <class Derived>
    struct is_solver_derived
    {
    private:
        template <int nx, int nv, int ni, int r>
        static std::true_type test(const Solver<Derived, nx, nv, ni, r> *);

        static std::false_type test(...);

    public:
        static constexpr bool value = decltype(test(std::declval<Derived *>()))::value;
    };

    template <class Derived>
    struct BFunc
    {
        static_assert(is_solver_derived<Derived>::value, "Derived must inherit from Solver");
        virtual void bf(Derived &solver, int inIx, int outIx) = 0;
        virtual void bfI(Derived &solver, int inIx, int outIx) = 0;
    };

    template <class Derived>
    struct BFuncOpen final : public BFunc<Derived>
    {
        static_assert(false, "Open boundary function not implemented for this solver");
    };

    template <class Derived>
    struct BFuncClosed final : public BFunc<Derived>
    {
        static_assert(false, "Closed boundary function not implemented for this solver");
    };

    template <class Derived>
    struct BFuncDriver final : public BFunc<Derived>
    {
        static_assert(false, "Driver boundary function not implemented for this solver");
    };
}