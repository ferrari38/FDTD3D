#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "pml3d.h"

typedef struct {
	int nx0, ny0, nz0, nx1, ny1, nz1;
	double ***exy, ***exz;
	double ***eyz, ***eyx;
	double ***ezx, ***ezy;
	double ***hxy, ***hxz;
	double ***hyz, ***hyx;
	double ***hzx, ***hzy;
	double ***aeyx, ***amyx;
	double ***aezx, ***amzx;
	double ***aexy, ***amxy;
	double ***aezy, ***amzy;
	double ***aexz, ***amxz;
	double ***aeyz, ***amyz;
	double ***beyx, ***bmyx;
	double ***bezx, ***bmzx;
	double ***bexy, ***bmxy;
	double ***bezy, ***bmzy;
	double ***bexz, ***bmxz;
	double ***beyz, ***bmyz;

} pml_var;

//PML吸収境界
int lpml = LPML; //PMLの次数
int order = ORDER; //PMLの層数
double rmax = RMAX; //要求精度[dB]

pml_var pml_x0, pml_x1, pml_y0, pml_y1, pml_z0, pml_z1;
double copml =  -1.5280063e-4;


//グローバル宣言
pml_var addpml(int nx0, int nx1, int ny0, int ny1, int nz0, int nz1);
pml_var epml(pml_var pml);
pml_var hpml(pml_var pml);
double ***allocate(int nx0, int nx1, int ny0, int ny1, int nz0, int nz1);
void afree(pml_var pml);


/*
 * @fn
 * PMLの初期設定
 */
void init_pml() {
	pml_x0 = addpml(0, LPML, 0, NY, 0, NZ);
	pml_x1 = addpml(NX-LPML, NX, 0, NY, 0, NZ);
	pml_y0 = addpml(0, NX, 0, LPML, 0, NZ);
	pml_y1 = addpml(0, NX, NY-LPML, NY, 0, NZ);
	pml_z0 = addpml(0, NX, 0, NY, 0, LPML);
	pml_z1 = addpml(0, NX, 0, NY, NZ-LPML, NZ);
}

/*
 * @fn
 * メモリの解放
 */
void free_pml() {
	afree(pml_x0);
	afree(pml_x1);
	afree(pml_y0);
	afree(pml_y1);
	afree(pml_z0);
	afree(pml_z1);
}

/*
 * @fn
 * メモリの解放
 */
void afree(pml_var pml) {
	free(pml.exy); free(pml.exz);
	free(pml.eyz); free(pml.eyx);
	free(pml.ezx); free(pml.ezy);

	free(pml.hxy); free(pml.hxz);
	free(pml.hyz); free(pml.hyx);
	free(pml.hzx); free(pml.hzy);

	free(pml.aeyx); free(pml.amyx);
	free(pml.aezx); free(pml.amzx);
	free(pml.aexy); free(pml.amxy);
	free(pml.aezy); free(pml.amzy);
	free(pml.aexz); free(pml.amxz);
	free(pml.aeyz); free(pml.amyz);

	free(pml.beyx); free(pml.bmyx);
	free(pml.bezx); free(pml.bmzx);
	free(pml.bexy); free(pml.bmxy);
	free(pml.bezy); free(pml.bmzy);
	free(pml.bexz); free(pml.bmxz);
	free(pml.beyz); free(pml.bmyz);
}

/*
 * @fn
 * 電界に対するPML
 */
 void e_pml() {
 	pml_x0 = epml(pml_x0);
 	pml_x1 = epml(pml_x1);
 	pml_y0 = epml(pml_y0);
 	pml_y1 = epml(pml_y1);
 	pml_z0 = epml(pml_z0);
 	pml_z1 = epml(pml_z1);
 }

 /*
 * @fn
 * 磁界に対するPML
 */
 void h_pml() {
 	pml_x0 = hpml(pml_x0);
 	pml_x1 = hpml(pml_x1);
 	pml_y0 = hpml(pml_y0);
 	pml_y1 = hpml(pml_y1);
 	pml_z0 = hpml(pml_z0);
 	pml_z1 = hpml(pml_z1);
 }


/*
 * @fn
 * 係数の計算
 */

pml_var addpml(int nx0, int nx1, int ny0, int ny1, int nz0, int nz1) {
	pml_var pml;

	pml.nx0 = nx0;
	pml.nx1 = nx1;
	pml.ny0 = ny0;
	pml.ny1 = ny1;
	pml.nz0 = nz0;
	pml.nz1 = nz1;

	pml.exy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.exz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.eyz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.eyx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.ezx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.ezy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.hxy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.hxz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.hyz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.hyx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.hzx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.hzy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);

	pml.aeyx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.aezx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.aexy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.aezy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.aexz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.aeyz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.amyx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.amzx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.amxy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.amzy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.amxz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.amyz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);

	pml.beyx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bezx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bexy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bezy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bexz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.beyz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bmyx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bmzx = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bmxy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bmzy = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bmxz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);
	pml.bmyz = allocate(nx0, nx1, ny0, ny1, nz0, nz1);

	for(int i=0; i<nx1; i++) {
		for(int j=0; j<ny1; j++) {
			for(int k=0; k<nz1; k++) {
				pml.exy[i][j][k] = 0.0;
				pml.exz[i][j][k] = 0.0;
				pml.eyz[i][j][k] = 0.0;
				pml.eyx[i][j][k] = 0.0;
				pml.ezx[i][j][k] = 0.0;
				pml.ezy[i][j][k] = 0.0;
				pml.hxy[i][j][k] = 0.0;
				pml.hxz[i][j][k] = 0.0;
				pml.hyz[i][j][k] = 0.0;
				pml.hyx[i][j][k] = 0.0;
				pml.hzx[i][j][k] = 0.0;
				pml.hzy[i][j][k] = 0.0;
			}
		}
	}
	

	double smax0x = copml*rmax*(order+1)/(lpml*DX);
	double smax0y = copml*rmax*(order+1)/(lpml*DY);
	double smax0z = copml*rmax*(order+1)/(lpml*DZ);

	double sigmxm, sigmxe, sigmym, sigmye, sigmzm, sigmze;

	for(int i=nx0; i<nx1-1; i++) {
		for(int j=ny0; j<ny1-1; j++) {
			for(int k=nz0; k<nz1-1; k++) {
				if(i<lpml) {
					sigmxm = pow((8-(double)i-0.5)/8,4)*smax0x;
					sigmxe = pow((8-(double)i)/8,4)*smax0x;
				}
				else if(i>=NX-lpml) {
					sigmxm = pow(((double)i-216+8+0.5)/8,4)*smax0x;
					sigmxe = pow(((double)i-216+8)/8,4)*smax0x;
				}
				else {
					sigmxm = 0.0;
					sigmxe = 0.0;
				}

				if(j<lpml) {
					sigmym = pow((8-(double)j-0.5)/8,4)*smax0y;
					sigmye = pow((8-(double)j)/8,4)*smax0y;
				}
				else if(j>=NY-lpml) {
					sigmym = pow(((double)j-216+8+0.5)/8,4)*smax0y;
					sigmye = pow(((double)j-216+8)/8,4)*smax0y;
				}
				else {
					sigmym = 0.0;
					sigmye = 0.0;
				}

				if(k<lpml) {
					sigmzm = pow((8-(double)k-0.5)/8,4)*smax0z;
					sigmze = pow((8-(double)k)/8,4)*smax0z;
				}
				else if(k>=NZ-lpml) {
					sigmzm = pow(((double)k-216+8+0.5)/8,4)*smax0z;
					sigmze = pow(((double)k-216+8)/8,4)*smax0z;
				}
				else {
					sigmzm = 0.0;
					sigmze = 0.0;
				}

				double epspml, mupml, sigmx, sigmy, sigmz, a;
				//Exの係数
				epspml = 0.25*(epsd[i][j][k]+epsd[i][j-1][k]+epsd[i][j][k-1]+epsd[i][j-1][k-1])*eps0;
				sigmy = sigmye*(epspml/eps0);
				sigmz = sigmze*(epspml/eps0);

				a = 0.5*sigmy*dt/epspml;
				pml.aexy[i][j][k] = (1.0-a)/(1.0+a);
				pml.bexy[i][j][k] = dt/epspml/(1.0+a)/dy;

				a = 0.5*sigmz*dt/epspml;
				pml.aexz[i][j][k] = (1.0-a)/(1.0+a);
				pml.bexz[i][j][k] = dt/epspml/(1.0+a)/dz;

				//Eyの係数
				epspml = 0.25*(epsd[i][j][k]+epsd[i-1][j][k]+epsd[i][j][k-1]+epsd[i-1][j][k-1])*eps0;
				sigmz = sigmze*(epspml/eps0);
				sigmx = sigmxe*(epspml/eps0);

				a = 0.5*sigmz*dt/epspml;
				pml.aeyz[i][j][k] = (1.0-a)/(1.0+a);
				pml.beyz[i][j][k] = dt/epspml/(1.0+a)/dz;

				a = 0.5*sigmx*dt/epspml;
				pml.aeyx[i][j][k] = (1.0-a)/(1.0+a);
				pml.beyx[i][j][k] = dt/epspml/(1.0+a)/dx;

				//Ezの係数
				epspml = 0.25*(epsd[i][j][k]+epsd[i-1][j][k]+epsd[i][j-1][k]+epsd[i-1][j-1][k])*eps0;
				sigmx = sigmxe*(epspml/eps0);
				sigmy = sigmye*(epspml/eps0);

				a = 0.5*sigmx*dt/epspml;
				pml.aezx[i][j][k] = (1.0-a)/(1.0+a);
				pml.bezx[i][j][k] = dt/epspml/(1.0+a)/dx;

				a = 0.5*sigmz*dt/epspml;
				pml.aezy[i][j][k] = (1.0-a)/(1.0+a);
				pml.bezy[i][j][k] = dt/epspml/(1.0+a)/dy;

				

				//Hxの係数
				mupml = 0.5*(mud[i][j][k]+mud[i-1][j][k])*mu0;
				epspml = 0.5*(epsd[i][j][k]+epsd[i-1][j][k])*eps0;
				sigmy = sigmym*(epspml/eps0);
				sigmz = sigmzm*(epspml/eps0);

				a = 0.5*sigmy*dt/epspml;
				pml.amxy[i][j][k] = (1.0-a)/(1.0+a);
				pml.bmxy[i][j][k] = dt/mupml/(1.0+a)/dy;

				a = 0.5*sigmz*dt/epspml;
				pml.amxz[i][j][k] = (1.0-a)/(1.0+a);
				pml.bmxz[i][j][k] = dt/mupml/(1.0+a)/dz;

				//Hyの係数
				mupml = 0.5*(mud[i][j][k]+mud[i][j-1][k])*mu0;
				epspml = 0.5*(epsd[i][j][k]+epsd[i][j-1][k])*eps0;
				sigmz = sigmzm*(epspml/eps0);
				sigmx = sigmxm*(epspml/eps0);

				a = 0.5*sigmz*dt/epspml;
				pml.amyz[i][j][k] = (1.0-a)/(1.0+a);
				pml.bmyz[i][j][k] = dt/mupml/(1.0+a)/dz;

				a = 0.5*sigmx*dt/epspml;
				pml.amyx[i][j][k] = (1.0-a)/(1.0+a);
				pml.bmyx[i][j][k] = dt/mupml/(1.0+a)/dx;

				//Hzの係数
				mupml = 0.5*(mud[i][j][k]+mud[i][j][k-1])*mu0;
				epspml = 0.5*(epsd[i][j][k]+epsd[i][j][k-1])*eps0;
				sigmx = sigmxm*(epspml/eps0);
				sigmy = sigmym*(epspml/eps0);

				a = 0.5*sigmx*dt/epspml;
				pml.amzx[i][j][k] = (1.0-a)/(1.0+a);
				pml.bmzx[i][j][k] = dt/mupml/(1.0+a)/dx;

				a = 0.5*sigmy*dt/epspml;
				pml.amzy[i][j][k] = (1.0-a)/(1.0+a);
				pml.bmzy[i][j][k] = dt/mupml/(1.0+a)/dy;

			}
		}
	}

	return pml;
}

/*
 * @fn
 * 電界の計算
 */
 pml_var epml(pml_var pml) {
	//Ex

	
	for(int i=pml.nx0; i<pml.nx1-1; i++) {
		for(int j=pml.ny0+1; j<pml.ny1-1; j++) {
			for(int k=pml.nz0+1; k<pml.nz1-1; k++) {
				pml.exy[i][j][k] = pml.aexy[i][j][k]*pml.exy[i][j][k]
					+pml.bexy[i][j][k]*(hz[i][j][k]-hz[i][j-1][k]);
				pml.exz[i][j][k] = pml.aexz[i][j][k]*pml.exz[i][j][k]
					+pml.bexz[i][j][k]*(hy[i][j][k-1]-hy[i][j][k]);
				ex[i][j][k] = pml.exy[i][j][k] + pml.exz[i][j][k];
			}
		}
	}
	



	//Ey
	for(int i=pml.nx0+1; i<pml.nx1-1; i++) {
		for(int j=pml.ny0; j<pml.ny1-1; j++) {
			for(int k=pml.nz0+1; k<pml.nz1-1; k++) {
				pml.eyz[i][j][k] = pml.aeyz[i][j][k]*pml.eyz[i][j][k]
					+pml.beyz[i][j][k]*(hx[i][j][k]-hx[i][j][k-1]);
				pml.eyx[i][j][k] = pml.aeyx[i][j][k]*pml.eyx[i][j][k]
					+pml.beyx[i][j][k]*(hz[i-1][j][k]-hz[i][j][k]);
				ey[i][j][k] = pml.eyz[i][j][k] + pml.eyx[i][j][k];
			}
		}
	}

	//Ez
	for(int i=pml.nx0+1; i<pml.nx1-1; i++) {
		for(int j=pml.ny0+1; j<pml.ny1-1; j++) {
			for(int k=pml.nz0; k<pml.nz1-1; k++) {
				pml.ezx[i][j][k] = pml.aezx[i][j][k]*pml.ezx[i][j][k]
					+pml.bezx[i][j][k]*(hy[i][j][k]-hy[i-1][j][k]);
				pml.ezy[i][j][k] = pml.aezy[i][j][k]*pml.ezy[i][j][k]
					+pml.bezy[i][j][k]*(hx[i][j-1][k]-hx[i][j][k]);
				ez[i][j][k] = pml.ezx[i][j][k] + pml.ezy[i][j][k];
			}
		}
	}

	return pml;

}

/*
 * @fn
 * 電界の計算
 */
pml_var hpml(pml_var pml) {
	//Hx
	for(int i=pml.nx0+1; i<pml.nx1-1; i++) {
		for(int j=pml.ny0; j<pml.ny1-1; j++) {
			for(int k=pml.nz0; k<pml.nz1-1; k++) {
				pml.hxy[i][j][k] = pml.amxy[i][j][k]*pml.hxy[i][j][k]
					+pml.bmxy[i][j][k]*(ez[i][j][k]-ez[i][j+1][k]);
				pml.hxz[i][j][k] = pml.amxz[i][j][k]*pml.hxz[i][j][k]
					+pml.bmxz[i][j][k]*(ey[i][j][k+1]-ey[i][j][k]);
				hx[i][j][k] = pml.hxy[i][j][k] + pml.hxz[i][j][k];
			}
		}
	}

	//Hy
	for(int i=pml.nx0; i<pml.nx1-1; i++) {
		for(int j=pml.ny0+1; j<pml.ny1-1; j++) {
			for(int k=pml.nz0; k<pml.nz1-1; k++) {
				pml.hyz[i][j][k] = pml.amyz[i][j][k]*pml.hyz[i][j][k]
					+pml.bmyz[i][j][k]*(ex[i][j][k]-ex[i][j][k+1]);
				pml.hyx[i][j][k] = pml.amyx[i][j][k]*pml.hyx[i][j][k]
					+pml.bmyx[i][j][k]*(ez[i+1][j][k]-ez[i][j][k]);
				hy[i][j][k] = pml.hyz[i][j][k] + pml.hyx[i][j][k];
			}
		}
	}

	//Hz
	for(int i=pml.nx0; i<pml.nx1-1; i++) {
		for(int j=pml.ny0; j<pml.ny1-1; j++) {
			for(int k=pml.nz0+1; k<pml.nz1-1; k++) {
				pml.hzx[i][j][k] = pml.amzx[i][j][k]*pml.hzx[i][j][k]
					+pml.bmzx[i][j][k]*(ey[i][j][k]-ey[i+1][j][k]);
				pml.hzy[i][j][k] = pml.amzy[i][j][k]*pml.hzy[i][j][k]
					+pml.bmzy[i][j][k]*(ex[i][j+1][k]-ex[i][j][k]);
				hz[i][j][k] = pml.hzx[i][j][k] + pml.hzy[i][j][k];
			}
		}
	}

	return pml;
}

/*
 * @fn
 * 動的配列の確保
 */
double ***allocate(int nx0, int nx1, int ny0, int ny1, int nz0, int nz1) {

	double ***a = malloc(nx1 * sizeof(double **));
	
	for(int i=0; i<nx1; i++) {
		a[i] = malloc(ny1 * sizeof(double*));
		for(int j=0; j<ny1; j++) {
			a[i][j] = malloc(nz1 * sizeof(double));
		}
	}

	return a;
}