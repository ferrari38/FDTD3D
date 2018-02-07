#include <stdio.h>
#include <math.h>
#include <complex.h>

#define K 5

double complex culcuEpsir(double ome);
double complex culcuX(int j, double ome);

double A[K+1];
double V[K+1];
double Omep;
double Omej[K+1];

double lam[50];
double complex epsir[50];

int main() {

	//データの代入
	A[0] = 0.760; A[1] = 0.024;
	A[2] = 0.010; A[3] = 0.071;
	A[4] = 0.601; A[5] = 4.384;

	V[0] = 0.053; V[1] = 0.241;
	V[2] = 0.345; V[3] = 0.870;
	V[4] = 2.494; V[5] = 2.214;

	Omep = 9.03;
	
	Omej[0] = 0.0; Omej[1] = 0.415;
	Omej[2] = 0.830; Omej[3] = 2.969;
	Omej[4] = 4.304; Omej[5] = 13.32;

	//データの横軸（λ）を作成
	for(int i=0; i<50; i++) {
		lam[i] = i*10 + 300;
	}

	//計算
	for(int i=0; i<50; i++) {
		double tmpOme = (4.135e-15*2.99e+8) / (lam[i]*1.0e-9);
		epsir[i] = culcuEpsir(tmpOme);
	}

	FILE *gp;
    gp = popen( "gnuplot -persist","w");
    fprintf(gp, "set yrange [%f:%f]\n", -30.0, 10.0);
    
    FILE *fpr = fopen("pointReal.dat", "w");
    for(int i=0; i<50; i++) {
    	fprintf(fpr, "%.1f %e\n", lam[i], creal(epsir[i]));
    }
    
    FILE *fpi = fopen("pointImag.dat", "w");
    for(int i=0; i<50; i++) {
    	fprintf(fpi, "%.1f %e\n", lam[i], cimag(epsir[i]));
    }
    
    fprintf(gp, "plot 'pointReal.dat' with lines\n");
    fprintf(gp, "replot 'pointImag.dat' with lines\n");
    
    fclose(fpr);
    fclose(fpi);
    
    fprintf(gp, "exit\n");
    pclose(gp);

	return 0;
}

/**
 * @fn
 * εを計算する関数
 *
 * @parm (ome) 式中ω
 * @return 計算結果
 */
 double complex culcuEpsir(double ome) {
	double complex ans = 0;
	double complex sigma = 0;

	for(int i=1; i<=K; i++) {
		sigma += culcuX(i, ome);
	}

	ans = 1 - A[0]*Omep*Omep / (ome*(ome+V[0]*I)) + sigma;

	return ans;
}

/**
 * @fn
 * xを計算する関数
 *
 * @parm (j) 式中添字j
 * @parm (ome) 式中ω
 * @return 計算結果
 */
double complex culcuX(int j, double ome) {
	double complex ans = 0;
	ans = A[j]*Omep*Omep / ((Omej[j]*Omej[j] - ome*ome) - ome*V[j]*I);

	return ans;
}