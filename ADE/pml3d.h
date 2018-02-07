void init_pml();
void free_pml();
void e_pml();
void h_pml();

extern double dt;
extern double dx; extern double dy; extern double dz;

extern double ex[NX][NY][NZ];
extern double ey[NX][NY][NZ];
extern double ez[NX][NY][NZ];

extern double hx[NX][NY][NZ];
extern double hy[NX][NY][NZ];
extern double hz[NX][NY][NZ];

extern double epsd[NX][NY][NZ];
extern double mud[NX][NY][NZ];

extern double eps0;
extern double mu0;