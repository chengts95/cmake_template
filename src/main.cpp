#include <stdio.h>
#include <RLC.h>
#include <matplotlib.hpp>
namespace plt = matplotlibcpp;
void swap(int & a, int & b)
{
    int temp;
    temp = b;
    b = a;
    a = temp;
}

int main(int argc, char const *argv[])
{
    RLCBranch rlc;
    std::vector<Real> data;
    for(int i=0;i<int(1e6);i++){
        data.emplace_back(sin(377.*i*1e-6));

    }
    plt::plot(data);
    plt::show();
    return 0;
}
