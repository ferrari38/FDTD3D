#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FDTD3D.h"

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
double compl =  -1.5290063e-4;

//グローバル宣言
void init_pml();
void addpml(pml_var pml, int nx0, int nx1, int ny0, int ny1, int nz0, int nz1);
void allocate(double ***a, int nx0, int nx1, int ny0, int ny1, int nz0, int nz1);

/*
 * @fn
 * メインプログラム
 *
 */
int main(int argc, const char * argv[]) {
	init_pml();
	return 0;
}

/*
 * @fn
 * PMLの初期設定
 */
void init_pml() {
	addpml(pml_x0, 0, LPML, 0, NY, 0, NZ);
	addpml(pml_x1, NX-LPML, NX, 0, NY, 0, NZ);
	/*
	addpml(pml_y0, 0, nx, 0, lpml, 0, nz);
	addpml(pml_y1, 0, nx, ny-lpml, ny, 0, nz);
	addpml(pml_z0, 0, nx, 0, ny, 0, lpml);
	addpml(pml_z1, 0, nx, 0, ny, nz-lpml, nz);
	*/
}

/*
 * @fn
 * 係数の計算
 */
void addpml(pml_var pml, int nx0, int nx1, int ny0, int ny1, int nz0, int nz1) {
	pml.nx0 = nx0;
	pml.nx1 = nx1;
	pml.ny0 = ny0;
	pml.ny1 = ny1;
	pml.nz0 = nz0;
	pml.nz1 = nz1;
	allocate(pml.exy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.exz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.eyz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.eyx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.ezx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.ezy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.hxy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.hxz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.hyz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.hyx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.hzx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.hzy, nx0, nx1, ny0, ny1, nz0, nz1);

	allocate(pml.aeyx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.aezx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.aexy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.aezy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.aexz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.aeyz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.amyx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.amzx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.amxy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.amzy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.amxz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.amyz, nx0, nx1, ny0, ny1, nz0, nz1);

	allocate(pml.beyx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bezx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bexy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bezy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bexz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.beyz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bmyx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bmzx, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bmxy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bmzy, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bmxz, nx0, nx1, ny0, ny1, nz0, nz1);
	allocate(pml.bmyz, nx0, nx1, ny0, ny1, nz0, nz1);

	for(int i=0; i<nx0; i++) {
		for(int j=0; j<ny0; j++) {
			for(int k=0; k<nz0; k++) {
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

	double smax0x = compl*rmax*(order+1)/(lpml+DX);
	double smax0y = compl*rmax*(order+1)/(lpml+DX);
	double smax0z = compl*rmax*(order+1)/(lpml+DX);

	double sigmxm, sigmxe, sigmym, sigmye, sigmzm, sigmze;

	for(int i=nx0; i<nx1-1; i++) {
		for(int j=ny0; j<ny1-1; j++) {
			for(int k=nz0; k<nz1-1; k++) {
				if(i<lpml) {
					sigmxm = pow((lpml-i-0.5)/lpml,order)*smax0x;
					sigmxe = pow((lpml-i)/lpml,order)*smax0x;
				}
				else if(i>=nx-lpml) {
					sigmxm = pow((i-nx+lpml+0.5)/lpml,order)*smax0x;
					sigmxe = pow((i-nx+lpml)/lpml,order)*smax0x;
				}
				else {
					sigmxm = 0.0;
					sigmxe = 0.0;
				}









			}
		}
	}

	double epsml = 0.25*(epsd)



}

/*
 * @fn
 * 動的配列の確保
 */
void allocate(double ***a, int nx0, int nx1, int ny0, int ny1, int nz0, int nz1) {
	a = malloc(nx1 * sizeof(double **));

	for(int i=0; i<nx1; i++) {
		a[i] = malloc(ny1 * sizeof(double*));
		for(int j=0; i<ny1; j++) {
			a[i][j] = malloc(nz1 * sizeof(double));
		}
	}
}