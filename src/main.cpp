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

template<int n,int lda>
int dgetrs(Real * a, Real *b) {

	
	for (int i = 0; i < n; i++) {

		int idx = IDX(i, 0, lda);
		for (int k = 0; k < i; k++)
			b[i] -= a[idx + k] * b[k];
	}

	for (int i = n - 1; i >= 0; i--) {
		int idx = IDX(i, 0, lda);
		for (int k = i + 1; k < n; k++)
			b[i] -= a[idx + k] * b[k];

		b[i] /= a[idx + i];
	}
	return 0;
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
