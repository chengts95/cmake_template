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
Real G[M * M] = {0};
Real b[M] = {0};
int main(int argc, char const *argv[])
{
    RLCBranch rlc;
    std::vector<Real> data;
    Real R1 = 5, R2 = 10;
    Real I2 = 1;
    Real U1 = 1.0, Rs = 1e-6;
    Real I1 = U1 / Rs;

    G[IDX(0, 0, M)] += 1.0 / R1;
    G[IDX(1, 2, M)] -= 1.0 / R1;
    G[IDX(2, 1, M)] -= 1.0 / R1;
    G[IDX(0, 0, M)] += 1.0 / Rs;
    G[IDX(1, 1, M)] += 1.0 / R2;
    G[IDX(1, 1, M)] += 1.0 / R1;
    b[0]+=I1;
    b[1]+=I2;
    return 0;
}
