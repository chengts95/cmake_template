#include <stdio.h>
#include <Components/RLC.h>

#include <Manager.h>
#include <Solver.h>
#include <omp.h>
#include <Fastor/Fastor.h>
#include <ctime>
#include <matplotlib.hpp>
//#define AoS
namespace plt = matplotlibcpp;
void swap(int &a, int &b)
{
    int temp;
    temp = b;
    b = a;
    a = temp;
}

Real A[3 * 3] = {0};
Fastor::Tensor<double, 2, 2> G;
using namespace Fastor;

template <int n, int lda>
inline void dgetrs(Real *a, Real *b)
{

    for (int i = 0; i < n; i++)
    {

        int idx = IDX(i, 0, lda);
        for (int k = 0; k < i; k++)
            b[i] -= a[idx + k] * b[k];
    }

    for (int i = n - 1; i >= 0; i--)
    {
        int idx = IDX(i, 0, lda);
        for (int k = i + 1; k < n; k++)
            b[i] -= a[idx + k] * b[k];

        b[i] /= a[idx + i];
    }
}
template <int n, int lda>
inline void dgetrf(Real *a)
{
    int min = n;
    for (int k = 0; k < min; k++)
    {
        for (int j = k + 1; j < n; j++)
        {

            a[IDX(j, k, lda)] /= a[IDX(k, k, lda)];
            for (int i = k + 1; i < n; i++)
            {

                a[IDX(j, i, lda)] -= a[IDX(j, k, lda)] * a[IDX(k, i, lda)];
            }
        }
    }
}

class UpdateCap
{
public:
    UpdateCap(int count, Real *sol, Real *b) : count(count), sol(sol), b(b)
    {
    }

    template <typename SOA_ACCESSOR>
    void operator()(SOA_ACCESSOR &Cap)
    {
        auto U = *sol;
        Real bIeq = 0;
        auto Geq = &Cap.Geq();
        auto Ieq = &Cap.Ieq();

//LibFlatArray::soa_accessor<Capacitor, 32, 32, 32, 0>::
//LibFlatArray::soa_accessor<Capacitor, 32, 32, 32, 0>::
#pragma omp simd reduction(+ \
                           : bIeq)
        for (int i = 0; i < count; i++)
        {
            	Real Ic = (U*Geq[i]+Ieq[i]);

            Ieq[i] = -(U * Geq[i] + Ic);
            bIeq += Ieq[i];
        }
        *b = -bIeq;
    }

private:
    int count;
    Real *sol;
    Real *b;
};

int main(int argc, char const *argv[])
{
    Resistor a{1, 0, 1};
    Resistor a1{1, 1, 2};
    Resistor a2{1, 0, 2};
    const int M = 1000;
    Capacitor c1[M];
    LibFlatArray::soa_vector<Capacitor> c2;

    std::vector<Real> data;

    Real b[2] = {0}, x[2] = {0};
    Tensor<Real, 2> B(b);
    Tensor<Real, 2> X(x);
    for (int i = 0; i < M; i++)
    {
        c1[i].C = 1e-3 / double(M);
        c1[i].port = {2, 0};
        c1[i].init_dt(1e-6);
        c2.push_back(c1[i]);
        c1[i].insert_MNA(A, 3);
    }
    Real Ieq = 0;
    UpdateCap capacitor(M, &B.data()[1], &Ieq);
    a.insert_MNA(A, 3);
    a1.insert_MNA(A, 3);

    //dgetrf<2, 3>(A + IDX(1, 1, 3));
    //dgetrf(0, 2, 2, A + IDX(1, 1, 3), 3);
    G(0, 0) = A[IDX(1, 1, 3)];
    G(0, 1) = A[IDX(1, 2, 3)];
    G(1, 0) = A[IDX(2, 1, 3)];
    G(1, 1) = A[IDX(2, 2, 3)];
    Tensor<double, 2, 2> C = inverse(G);

    auto t_start = std::chrono::high_resolution_clock::now();
//double t0 = omp_get_wtime();
#ifdef AoS
    for (int i = 0; i < int(1e6); i++)
    {
        //c2.callback(capacitor);
        Ieq = 0;
#pragma omp simd reduction(+ \
                           : Ieq)
        for (int i = 0; i < M; i++)
        {
            c1[i].compute_Ieq(B(1));
            Ieq += c1[i].Ieq;
        }
        B.data()[0] = 100;
        B.data()[1] = -Ieq;
        //dgetrs<2, 2>(G.data(),B.data());
        //B = solve<SolveCompType::SimpleInv>(G, B);
        B = matmul(C, B);

        data.emplace_back(B(1));
    }
#else
    for (int i = 0; i < int(1e6); i++)
    {
        c2.callback(capacitor);
        B.data()[0] = 100;
        B.data()[1] = Ieq;
        //dgetrs<2, 2>(G.data(),B.data());
        //B = solve<SolveCompType::SimpleInv>(G, B);
        B = matmul(C, B);

        data.emplace_back(B(1));
    }
#endif
    auto t_end = std::chrono::high_resolution_clock::now();
    printf("SoA: %f ms\n", std::chrono::duration<double, std::milli>(t_end - t_start).count());
    plt::plot(data);
    plt::show();
    return 0;
}
