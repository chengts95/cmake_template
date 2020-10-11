#pragma once
#include "common.h"
#include "Components/IComponent.h"
#include <array>
using int2 = std::array<int,2>;

struct Resistor 
{
	Real R = 0;
	int2 port;
	template <typename T>
	void insert_MNA(T &A, int N);
	
};

struct Capacitor
{
	Real C = 0;
	int2 port;
	Real Ieq=0,Geq;
	void init_dt(Real dt);
	template <typename T>
	void insert_MNA(T &A, int N);
	void compute_Ieq(Real U);
};
struct Inductor
{
	Real L = 0;
	int2 port;
	Real Ieq=0,Geq;
	void init_dt(Real dt);
	template <typename T>
	void insert_MNA(T &A, int N);
	void compute_Ieq(Real U);
};


template <typename T>
inline void Resistor::insert_MNA(T &A, int N)
{
	insert_Y(A, N, port[0], port[1], 1. / R);
}

template <typename T>
inline void Capacitor::insert_MNA(T &A, int N) {
    insert_Y(A, N, port[0], port[1], Geq);
}

inline void Capacitor::compute_Ieq(Real U) {
	Real Ic = (U*Geq+Ieq);
    Ieq = -(U*Geq+Ic);
}
inline void Capacitor::init_dt(Real dt) {
    Geq = (2.0 * C) / dt;
}


inline void Inductor::init_dt(Real dt) {
    Geq = dt / (2.0 * L);
}

template <typename T>
inline void Inductor::insert_MNA(T &A, int N) {
    insert_Y(A, N, port[0], port[1], Geq);
}

inline void Inductor::compute_Ieq(Real U) {
	Real Il = U*Geq+Ieq;
    Ieq = U*Geq+Il;
}

