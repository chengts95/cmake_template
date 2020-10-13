#include <stdio.h>
#include <RLC.h>
#include <matplotlib.hpp>
namespace plt = matplotlibcpp;
void swap(int &a, int &b)
{
    int temp;
    temp = b;
    b = a;
    a = temp;
}
const int M = 2;
#define IDX(i, j, lda) (i * lda + j)

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

Real G[M * M] = {0};
Real b[M] = {0};
int main(int argc, char const *argv[])
{
    const Real dt =1e-6;
    RLCBranch rlc;
    std::vector<Real> data;
    Real R1 = 500, R2 = 10;
    Real C2 = 1e-3;
    Real Geq = 2.0*C2/dt;
    Real Ieq = 0;
    Real I2 = 1;
    Real U1 = 1.0, Rs = 1e-6;
    Real I1 = U1 / Rs;
    std::vector<Real> u;
    G[IDX(0, 0, M)] += 1.0 / R1;
    G[IDX(1, 0, M)] -= 1.0 / R1;
    G[IDX(0, 1, M)] -= 1.0 / R1;
    G[IDX(0, 0, M)] += 1.0 / Rs;
    G[IDX(1, 1, M)] +=  Geq;
    G[IDX(1, 1, M)] += 1.0 / R1;
    dgetrf<2,2>(G);
    for(int i=0;i<1000000;i++){
        b[0]=I1;
        b[1]=-Ieq;
        dgetrs<2,2>(G,b);
        Real U = b[1];
        Real Ic = U*Geq+Ieq;
        Ieq =-Geq-Ic;
        u.emplace_back(U);
    }
 
    plt::plot(u);
    plt::show();
    printf("%f %f\n",b[0],b[1]);
    return 0;
}
