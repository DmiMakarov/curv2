#define pi  3.14159265359

typedef struct{
    double* w1_;
    double* w1;
    double* w2_;
    double* w2;
    double* curv;
    double* dualCurv;
    double* trueCurv;
    double* lens;
}curvs;

void free_curvs(curvs* curv);
double my_abs(double x);
void printp(double * point);
double len(double* p1, double* p2);

void get_point_by_index(double* point,double* points, int indx);
double* get_dual_points(int* faces, double *points, double* normals,int numFaces);
int* build_dual_faces(double* points,double* normals,double* dualPoints, int numPoints, int numDualPoints);
int find_dual_face(double* point, double* dualPoints, int* dualFaces, int numDualPoints, int numDualFaces);
int find_nearest_face(double* point, double* points, int* faces,int numPoints,int numFaces);

void circle(double* cir,double* param);
void Dcircle(double* p,double *param);
void DDcircle(double* p,double *param);
double circle_curv(double *param);

void spiral(double* cir,double* param);
void Dspiral(double* p,double *param);
void DDspiral(double* p,double *param);
double spiral_curv(double *param);

void ellipse(double* cir,double* param);
void Dellipse(double* p,double *param);
void DDellipse(double* p,double *param);
double ellipse_curv(double *param);

void delta(double* cir,double* param);
void Ddelta(double* p,double *param);
void DDdelta(double* p,double *param);
double delta_curv(double *param);

void loop(double* cir,double* param);
void Dloop(double* p,double *param);
void DDloop(double* p,double *param);
double loop_curv(double *param);

void parabola(double* cir,double* param);
void Dparabola(double* p,double *param);
void DDparabola(double* p,double *param);
double parabola_curv(double *param);

curvs* compute_curvs(double *points, double* normals, int *faces, double *dualPoints,int *dualFaces, int numPoints, int numFaces,void (*fig)(double*,double *), void (*Dfig)(double*,double*),void (*DDfig)(double*,double*), double (*figCurv)(double*), double* param,int numFigPoints, double t0, double tmax,int IsFlat);