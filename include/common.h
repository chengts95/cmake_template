#pragma once
#include <vector>

#define ROW_IDX(i, j, ld) (((i) * (ld)) + (j))
#define COL_IDX(i, j, ld) (((j) * (ld)) + (i))
#define IDX(i, j, ld) ROW_IDX(i, j, ld)
typedef double Real;
template<class T>
struct Triplet{
	int i;
	int j;
	T v;
};

typedef std::vector<Triplet<Real>> COO;

inline void insert_Y(Real *MNA, int dim, int n1, int n2, Real geq)
{

	MNA[IDX(n2, n2, dim)] += geq;
	MNA[IDX(n1, n1, dim)] += geq;
	MNA[IDX(n1, n2, dim)] -= geq;
	MNA[IDX(n2, n1, dim)] -= geq;
}

inline void insert_Y(COO &MNA, int dim, int n1, int n2, Real geq)
{
	if (n1 > 0)
	{
		MNA.emplace_back(Triplet<Real>{n1, n1, geq});
		if (n2 > 0)
		{
			MNA.emplace_back(Triplet<Real>{n1, n2, -geq});
			MNA.emplace_back(Triplet<Real>{n2, n1, -geq});
			MNA.emplace_back(Triplet<Real>{n2, n2, geq});
		}
	}
	else if (n2 > 0)
	{
		MNA.emplace_back(Triplet<Real>{n2, n2, geq});
	}
}