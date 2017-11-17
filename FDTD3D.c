#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "FDTD3D.h"

int nstep = NSTEP; //計算ステップ数
double dx = DX; double dy = DY; double dz = DZ;
double dt, t; //時間ステップ、時間

//電界、磁界の配列、計算の配列
double ex[NX][NY][NZ], ey[NX][NY][NZ], ez[NX][NY][NZ];
double hx[NX][NY][NZ], hy[NX][NY][NZ], hz[NX][NY][NZ];
double aex[NX][NY][NZ], aey[NX][NY][NZ], aez[NX][NY][NZ]; //係数
double bexy[NX][NY][NZ], bexz[NX][NY][NZ]; //係数
double beyx[NX][NY][NZ], beyz[NX][NY][NZ]; //係数
double bezx[NX][NY][NZ], bezy[NX][NY][NZ]; //係数
double amx[NX][NY][NZ], amy[NX][NY][NZ], amz[NX][NY][NZ]; //係数
double bmxy[NX][NY][NZ], bmxz[NX][NY][NZ]; //係数
double bmyx[NX][NY][NZ], bmyz[NX][NY][NZ]; //係数
double bmzx[NX][NY][NZ], bmzy[NX][NY][NZ]; //係数

//媒質定数の配列と背景媒質定数
double epsd[NX][NY][NZ]; double sgmed[NX][NY][NZ];
double mud[NX][NY][NZ]; double sgmmd[NX][NY][NZ];
double epsbk = 1.0;
double mubk = 1.0;
double sigebk = 0.0;
double sigmbk = 0.0;

//散乱体
int ic = NX/2; //散乱直方体の中心
int jc = NY/2;
int kc = NZ/2;
int lx2 = 20; //直方体の寸法/2
int ly2 = 20;
int lz2 = 20;
double epsr = 10.0; //直方体の比誘電率

//励振パルス
double duration; //パルス幅
double t0; //ピーク時刻
double dl = 0.001; //ハード給電の間隙
int ifed, jfed, kfed;

//定数
double c = 2.9979246e+8; //光速
double eps0 = 8.854188e-12; //真空の誘電率
double mu0 = 1.256637e-6; //真空の透磁率


//グローバル宣言
void setup();
void epsmu();
void feed();
void e_cal();
void h_cal();

/*
 * @fn
 * メインプログラム
 *
 */
int main(int argc, const char * argv[]) {

    setup(); //FDTDの初期設定

    FILE *fp = fopen("eztm_epsr.dat", "w");

    t = dt;
    for(int n=1; n<nstep; n++) { //繰り返し計算
        printf("Time step: %d\n", n);
        e_cal(); //電界の計算
        //電界に対するPML
        feed();//電流源の励振
        t = t + 0.5*dt; //時間の更新
        h_cal(); //磁界の計算
        //磁界に対するPML
        t = t + 0.5*dt; //時間の更新
        //計算結果の出力
        fprintf(fp,"%e %e\n", t, ez[70][100][100]);
    }

    fclose(fp);

    return 0;
}


/*
 * @fn
 * FDTDの初期設定をする関数
 */
 void setup() {
    double eps, sgm, mu, sgmm, a;

    //時間ステップ
    dt = 0.99999/(c*sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));

    //背景媒質
    for(int i=0; i<NX; i++) {
        for(int j=0; j<NY; j++) {
            for(int k=0; k<NZ; k++) {
                epsd[i][j][k] = epsbk;
                mud[i][j][k] = mubk;
                sgmed[i][j][k] = sigebk;
                sgmmd[i][j][k] = sigmbk;
            }
        }
    }
    //散乱直方体の媒質定数（背景媒質に上書き）
    epsmu();

    //波元の位置における係数
    /////////未完成

    //係数の計算
    for(int i=0; i<NX; i++) {
        for(int j=0; j<NY; j++) {
            for(int k=0; k<NZ; k++) {

                eps = 0.25*(epsd[i][j][k]+epsd[i][j-1][k]+epsd[i][j][k-1]+epsd[i][j-1][k-1])*eps0;
                sgm = 0.25*(sgmed[i][j][k]+sgmed[i][j-1][k]+sgmed[i][j][k-1]+sgmed[i][j-1][k-1]);
                a = 0.5*sgm*dt/eps;
                aex[i][j][k] = (1.0-a)/(1.0+a);
                bexy[i][j][k] = dt/eps/(1.0+a)/dy;
                bexz[i][j][k] = dt/eps/(1.0+a)/dz;

                eps = 0.25*(epsd[i][j][k]+epsd[i-1][j][k]+epsd[i][j][k-1]+epsd[i-1][j][k-1])*eps0;
                sgm = 0.25*(sgmed[i][j][k]+sgmed[i-1][j][k]+sgmed[i][j][k-1]+sgmed[i-1][j][k-1]);
                a = 0.5*sgm*dt/eps;
                aey[i][j][k] = (1.0-a)/(1.0+a);
                beyx[i][j][k] = dt/eps/(1.0+a)/dx;
                beyz[i][j][k] = dt/eps/(1.0+a)/dz;

                eps = 0.25*(epsd[i][j][k]+epsd[i-1][j][k]+epsd[i][j-1][k]+epsd[i-1][j-1][k])*eps0;
                sgm = 0.25*(sgmed[i][j][k]+sgmed[i-1][j][k]+sgmed[i][j-1][k]+sgmed[i-1][j-1][k]);
                a = 0.5*sgm*dt/eps;
                aez[i][j][k] = (1.0-a)/(1.0+a);
                bezx[i][j][k] = dt/eps/(1.0+a)/dx;
                bezy[i][j][k] = dt/eps/(1.0+a)/dy;

                mu = 0.5*(mud[i][j][k]+mud[i-1][j][k])*mu0;
                sgmm = 0.5*(sgmmd[i][j][k]+sgmmd[i-1][j][k]);
                a = 0.5*sgmm*dt/mu;
                amx[i][j][k] = (1.0-a)/(1.0+a);
                bmxy[i][j][k] = dt/mu/(1.0+a)/dy;
                bmxz[i][j][k] = dt/mu/(1.0+a)/dz;

                mu = 0.5*(mud[i][j][k]+mud[i][j-1][k])*mu0;
                sgmm = 0.5*(sgmmd[i][j][k]+sgmmd[i][j-1][k]);
                a = 0.5*sgmm*dt/mu;
                amy[i][j][k] = (1.0-a)/(1.0+a);
                bmyz[i][j][k] = dt/mu/(1.0+a)/dz;
                bmyx[i][j][k] = dt/mu/(1.0+a)/dx;

                mu = 0.5*(mud[i][j][k]+mud[i][j][k-1])*mu0;
                sgmm = 0.5*(sgmmd[i][j][k]+sgmmd[i][j][k-1]);
                a = 0.5*sgmm*dt/mu;
                amz[i][j][k] = (1.0-a)/(1.0+a);
                bmzx[i][j][k] = dt/mu/(1.0+a)/dx;
                bmzy[i][j][k] = dt/mu/(1.0+a)/dy;
            }
        }
    }

    //電界磁界の初期化
    for(int i=0; i<NX; i++) {
        for(int j=0; j<NY; j++) {
            for(int k=0; k<NZ; k++) {
                ex[i][j][k] = 0.0;
                ey[i][j][k] = 0.0;
                ez[i][j][k] = 0.0;

                hx[i][j][k] = 0.0;
                hy[i][j][k] = 0.0;
                hz[i][j][k] = 0.0;
            }
        }
    }
 }

 /*
  * @fn
  * 直方体の媒質定数
  */
void epsmu() {
    for(int i=ic-lx2; i<ic+lx2; i++) {
        for(int j=jc-ly2; j<jc+ly2; j++) {
            for(int k=kc-lz2; k<kc+lz2; k++) {
                epsd[i][j][k] = epsr;
                mud[i][j][k] = 1.0;
                sgmed[i][j][k] = 0.0;
                sgmmd[i][j][k] = 0.0;
            }
        }
    }
}

 /*
  * @fn
  * 励振波源
  */
void feed() {
    ifed = 60;
    jfed = 100;
    kfed = 100;
    duration = 0.1e-9; //パルス幅
    t0 = 4.0 * duration; //ピーク時刻
    ez[ifed][jfed][kfed] = exp(-pow((t-t0)/duration,2.0))/dl; //ハード給電
}


 /*
  * @fn
  * 電界の計算をする関数
  */
void e_cal() {
    //Ex
    for(int i=0; i<NX-1; i++) {
        for(int j=1; j<NY-1; j++) {
            for(int k=1; k<NZ-1; k++) {
                ex[i][j][k]
                = aex[i][j][k]*ex[i][j][k]
                + bexy[i][j][k]*(hz[i][j][k]-hz[i][j-1][k])
                - bexz[i][j][k]*(hy[i][j][k]-hy[i][j][k-1]);
            }
        }
    }

    //Ey
    for(int i=1; i<NX-1; i++) {
        for(int j=0; j<NY-1; j++) {
            for(int k=1; k<NZ-1; k++) {
                ey[i][j][k]
                = aey[i][j][k]*ey[i][j][k]
                + beyz[i][j][k]*(hx[i][j][k]-hx[i][j][k-1])
                - beyx[i][j][k]*(hz[i][j][k]-hz[i-1][j][k]);
            }
        }
    }

    //Ez
    for(int i=1; i<NX-1; i++) {
        for(int j=1; j<NY-1; j++) {
            for(int k=0; k<NZ-1; k++) {
                ez[i][j][k]
                = aez[i][j][k]*ez[i][j][k]
                + bezx[i][j][k]*(hy[i][j][k]-hy[i-1][j][k])
                - bezy[i][j][k]*(hx[i][j][k]-hx[i][j-1][k]);
            }
        }
    }
}

/*
 * @fn
 * 磁界の計算
 *
 */
 void h_cal() {
    //Hx
    for(int i=1; i<NX-1; i++) {
        for(int j=0; j<NY-1; j++) {
            for(int k=0; k<NZ-1; k++) {
                hx[i][j][k]
                = amx[i][j][k]*hx[i][j][k]
                - bmxy[i][j][k]*(ez[i][j+1][k]-ez[i][j][k])
                + bmxz[i][j][k]*(ey[i][j][k+1]-ey[i][j][k]);
            }
        }
    }

    //Hy
    for(int i=0; i<NX-1; i++) {
        for(int j=1; j<NY-1; j++) {
            for(int k=0; k<NZ-1; k++) {
                hy[i][j][k]
                = amy[i][j][k]*hy[i][j][k]
                - bmyz[i][j][k]*(ex[i][j][k+1]-ex[i][j][k])
                + bmyx[i][j][k]*(ez[i+1][j][k]-ez[i][j][k]);
            }
        }
    }

    //Hz
    for(int i=0; i<NX-1; i++) {
        for(int j=0; j<NY-1; j++) {
            for(int k=1; k<NZ-1; k++) {
                hz[i][j][k]
                = amz[i][j][k]*hz[i][j][k]
                - bmzx[i][j][k]*(ey[i+1][j][k]-ey[i][j][k])
                + bmzy[i][j][k]*(ex[i][j+1][k]-ex[i][j][k]);
            }
        }
    }
}