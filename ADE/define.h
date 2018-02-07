//解析領域
#define NXX 200 //解析領域分割数
#define NXY 200 //解析領域分割数
#define NXZ 200 //解析領域分割数

#define NSTEP 300 //計算ステップ数
#define DX 0.005
#define DY 0.005
#define DZ 0.005

#define K 5

//PML吸収境界
#define LPML 8 //PMLの次数
#define ORDER 4 //PMLの層数
#define RMAX -120.0 //要求精度[dB]

//全計算領域
#define NX (NXX + 2*LPML)
#define NY (NXY + 2*LPML)
#define NZ (NXZ + 2*LPML)