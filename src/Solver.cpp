
#include "Solver.h"
#include "string.h"



#define MIN(x,y)  (((x) < (y)) ? (x) : (y))

#define fabs(x) ((x)<0?-(x):(x))

inline void swap(Real & a, Real & b) {

	register Real temp;
	temp = a;
	a = b;
	b = temp;
	

};

inline void swapline(Real *a , Real * b,int size) {

	Real * temp=new Real[size];
	memcpy(temp,a,sizeof(Real)*size);
	memcpy(a, b, sizeof(Real)*size);
	memcpy(b, temp, sizeof(Real)*size);
	delete [] temp;
};

inline void dlaswp(Real *a, int m,int n,int lda,int * ipiv) {

	Real * temp = new Real[n];
	int rowsize = sizeof(Real)*n;
	int ip ;
	for (int i = 0; i < m;i++) {
		ip = ipiv[i];
		if(ip!=i){
			Real * idx = a+IDX(i, 0, lda);
			Real * idx2 = a+IDX(ip, 0, lda);
			memcpy(temp, idx, sizeof(Real)*n);
			memcpy(idx, idx2, sizeof(Real)*n);
			memcpy(idx2, temp, sizeof(Real)*n);
		}
	
	}

	delete[] temp;
};

int dgetrf(int matrix_layout, int m, int n, Real * a, int lda)
{
	int min = MIN(m, n);
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

	return 0;
}

int dgetrf_pivot(int matrix_layout, int m, int n, Real * a, int lda,int * ipiv)
{
	
	Real Amax, Acur;
	int pcur;
	int j;
	int min = MIN(m, n);
	for (int i = 0; i < m; i++)
		ipiv[i] = i; //Unit permutation matrix, P[N] initialized with N

	for (int k = 0; k < min; k++) {
		Amax = 0.0;
		pcur = k;
		for (int i = k; i < m; i++) {
			Acur = fabs(a[IDX(i, k, lda)]);
			if (Acur > Amax) {
				Amax = Acur;
				pcur = i;
			}
		
		}
		if (pcur != k) {
			
			j = ipiv[k];
			ipiv[k] = ipiv[pcur];
			ipiv[pcur] = j;

			swapline(a+IDX(pcur,0,lda), a + IDX(j, 0, lda),n);
			
		}

		for (j = k + 1; j < m; j++) {
			a[IDX(j,k,lda)] /= a[IDX(k, k, lda)];
		}

		for (j = k + 1; j < m; j++){
			for (int i = k + 1; i < n; i++){
				a[IDX(j, i, lda)] -= a[IDX(j, k, lda)] * a[IDX(k, i, lda)];
		}
	}
	}


	return 0;
}

int dgetrs(int matrix_layout, int n, int rhs, Real * a, int lda,  Real *b, int ldb) {

	
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

int dgetrs_pivot(int matrix_layout, int n, int rhs, Real * a, int lda, int * ipiv,Real *b,int ldb) {
	
	dlaswp(b,n,rhs,ldb,ipiv);
	for (int i = 0; i < n; i++) {
		
		int idx = IDX(i, 0, lda);
		for (int k = 0; k < i; k++)
			b[i] -= a[idx+k] * b[k];
	}

	for (int i = n - 1; i >= 0; i--) {
		int idx = IDX(i, 0, lda);
		for (int k = i + 1; k < n; k++)
			b[i] -= a[idx+k] * b[k];

		b[i]  /= a[idx+i];
	}
	return 0;
}
int dgetrf_pivot2(int matrix_layout, int m, int n, Real * a, int lda, int * ipiv)
{

	Real Amax, Acur;
	int pcur;

	for (int i = 0; i < m; i++)
		ipiv[i] = i; //Unit permutation matrix, P[N] initialized with N

	for (int j = 0; j < n; j++) {
		Amax = -1, pcur = j;

		for (int i = 0; i < m; i++) {
			bool temp = i < j;
			int idx = IDX(i, 0, lda);
			if (temp) {
				
				swap(a[IDX(i, j, lda)], a[IDX(ipiv[i], j, lda)]);
			}
		
			int min = MIN(j, i)-1;
			
			for(int k=0;k<min;k++){
				
				a[idx+j] = a[idx+j] - a[idx + k] * a[IDX(k, j, lda)];
			}

			if (temp) {
				a[idx+j] /= a[idx+i];
			}else {
				
				if ((Acur = fabs(a[idx + j])) >Amax) {
					Amax = Acur;
					pcur = i;
				
				}
			}

		}
		if (Amax >= 0) {
		
			ipiv[j] = pcur;
			if (pcur != j) {
				
				for (int k = 0; k < j; k++) {
					
					swap(a[IDX(j,k, lda)], a[IDX(pcur,k, lda)]);
				}
			}
		
		}
	}


	return 0;
}

//void LUPInvert(Real *A, int *P, int N, Real *IA) {
//
//	for (int j = 0; j < N; j++) {
//		for (int i = 0; i < N; i++) {
//			if (P[i] == j)
//				IA[i][j] = 1.0;
//			else
//				IA[i][j] = 0.0;
//
//			for (int k = 0; k < i; k++)
//				IA[i][j] -= A[i][k] * IA[k][j];
//		}
//
//		for (int i = N - 1; i >= 0; i--) {
//			for (int k = i + 1; k < N; k++)
//				IA[i][j] -= A[i][k] * IA[k][j];
//
//			IA[i][j] = IA[i][j] / A[i][i];
//		}
//	}
//}