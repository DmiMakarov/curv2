#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "curvature.h"

//param[0] = t
//param[1] = r
void circle(double*cir, double* param){
    cir[2] = 0;
    cir[0] = param[1]*sin(param[0]);
    cir[1] = param[1]*cos(param[0]);
}

void Dcircle(double*p,double *param){
    p[0] = param[1]*cos(param[0]);
    p[1] = -param[1]*sin(param[0]);
    p[2] = 0;
}

void DDcircle(double*p,double *param){
    p[0] = -param[1]*sin(param[0]);
    p[1] = -param[1]*cos(param[0]);
    p[2] = 0;
}

double circle_curv(double *param){
    return 1/param[1];
}

//param[0] = t
//param[1] = k 
void spiral(double*cir, double* param){
    cir[2] = 0;
    cir[0] = param[1]*param[0]*sin(param[0]);
    cir[1] = param[1]*param[0]*cos(param[0]);
}

void Dspiral(double*p,double *param){
    p[0] = param[1]*sin(param[0]) + param[1]*param[0]*cos(param[0]);
    p[1] = param[1]*cos(param[0]) - param[1]*param[0]*sin(param[0]);
    p[2] = 0;
}

void DDspiral(double*p,double *param){
    p[0] = 2*param[1]*cos(param[0]) - param[1]*param[0]*sin(param[0]);
    p[1] = -2*param[1]*sin(param[0]) - param[1]*param[0]*cos(param[0]);
    p[2] = 0;
}

double spiral_curv(double *param){
    //(phi**2 + 2)/(np.sqrt((phi**2 + 1)**3)*k)
    double curv = (param[0]*param[0] + 2);
    curv /= sqrt(param[0]*param[0]+1);
    curv /= sqrt(param[0]*param[0]+1);
    curv /= sqrt(param[0]*param[0]+1);
    curv /= param[1];
    return curv;
}

//param[0] = t
//param[1] = a
//param[2] = b 
void ellipse(double*cir, double* param){
    cir[2] = 0;
    cir[0] = param[1]*sin(param[0]);
    cir[1] = param[2]*cos(param[0]);
}

void Dellipse(double*p,double *param){
    p[0] = param[1]*cos(param[0]);
    p[1] = -param[2]*sin(param[0]);
    p[2] = 0;
}

void DDellipse(double*p,double *param){
    p[0] = -param[1]*sin(param[0]);
    p[1] = -param[2]*cos(param[0]);
    p[2] = 0;
}

double ellipse_curv(double *param){
    //a*b/np.sqrt((a**2*(np.cos(phi[:-1]))**2 + b**2*(np.sin(phi[:-1]))**2)**3)
    double curv = param[1]*param[2];
    curv /= sqrt(param[1]*param[1]*cos(param[0])*cos(param[0]) + param[2]*param[2]*sin(param[0])*sin(param[0]));
    curv /= sqrt(param[1]*param[1]*cos(param[0])*cos(param[0]) + param[2]*param[2]*sin(param[0])*sin(param[0]));
    curv /= sqrt(param[1]*param[1]*cos(param[0])*cos(param[0]) + param[2]*param[2]*sin(param[0])*sin(param[0]));
    return curv;
}

//param[0] = t
//param[1] = a
//param[2] = b
void delta(double*cir, double* param){
    if(my_abs(1 - cos(param[2]*param[0]/param[1]))<1e-5){
        param[0] +=0.0001;
    }
    cir[0] = (param[2] - param[1])*cos(param[0]) + param[1]*cos((param[2] - param[1])*param[0]/param[1]);
    cir[1] = (param[2] - param[1])*sin(param[0]) - param[1]*sin((param[2] - param[1])*param[0]/param[1]);
    cir[2] = 0;
}

void Ddelta(double*p,double *param){
    if(my_abs(1 - cos(param[2]*param[0]/param[1]))<1e-5){
        param[0] +=0.0001;
    }
    p[0] = -(param[2] - param[1])*sin(param[0]) - (param[2] - param[1])*sin((param[2] - param[1])*param[0]/param[1]);
    p[1] = (param[2] - param[1])*cos(param[0]) - (param[2] - param[1])*cos((param[2] - param[1])*param[0]/param[1]);
    p[2] = 0;
}

void DDdelta(double*p,double *param){
    if(my_abs(1 - cos(param[2]*param[0]/param[1]))<1e-5){
        param[0] +=0.0001;
    }
    p[0] = -(param[2] - param[1])*cos(param[0]) - (param[2] - param[1])*(param[2] - param[1])*cos((param[2] - param[1])*param[0]/param[1])/param[1];
    p[1] = -(param[2] - param[1])*sin(param[0]) + (param[2] - param[1])*(param[2] - param[1])*sin((param[2] - param[1])*param[0]/param[1])/param[1];
    p[2] = 0;
}

double delta_curv(double *param){
    //((b-a)**2-((b-a)**3)/a) * (1 - np.cos(b*phi[:-1]/a))/(np.sqrt(2*(b-a)**2)*(1-np.cos(b*phi[:-1]/a)))**3
    if(my_abs(1 - cos(param[2]*param[0]/param[1]))<1e-5){
        param[0] +=0.0001;
    }
    double curv = (param[2] - param[1])*(param[2] - param[1])*(1 - (param[2] - param[1])/param[1]);
    curv *=(1 - cos(param[2]*param[0]/param[1]));
    curv /= sqrt(2*(param[2] - param[1])*(param[2] - param[1]))*(1 - cos(param[2]*param[0]/param[1]));
    curv /= sqrt(2*(param[2] - param[1])*(param[2] - param[1]))*(1 - cos(param[2]*param[0]/param[1]));
    curv /= sqrt(2*(param[2] - param[1])*(param[2] - param[1]))*(1 - cos(param[2]*param[0]/param[1]));
    return curv;
}

//param[0] = t
void loop(double*cir, double* param){
    if(my_abs(param[0])<1e-5 || my_abs(param[0]-2)<1e-5){
        param[0] +=0.0001;
    }
    cir[0] = param[0]*(2-param[0]);
    cir[1] = param[0]*param[0]*(2-param[0]);
    cir[2] = 0;
}

void Dloop(double*p,double *param){
    if(my_abs(param[0])<1e-5 || my_abs(param[0]-2)<1e-5){
        param[0] +=0.0001;
    }
    p[0] = 2*(1 - param[0]);
    p[1] = 4*param[0] - 3*param[0]*param[0];
    p[2] = 0;
}

void DDloop(double*p,double *param){
    if(my_abs(param[0])<1e-5 || my_abs(param[0]-2)<1e-5){
        param[0] +=0.0001;
    }
    p[0] = -2;
    p[1] = 4 - 6*param[0];
    p[2] = 0;
}

double loop_curv(double *param){
    // 2*(3*t**3 - 6*t+4)/(np.sqrt(9*t**4 -24*t**3 + 20*t**2 -8*t +4))**3
    if(my_abs(param[0])<1e-5 || my_abs(param[0]-2)<1e-5){
        param[0] +=0.0001;
    }
    double curv = 2*(3*param[0]*param[0]*param[0]-6*param[0]+4);
    curv /= sqrt(9*param[0]*param[0]*param[0]*param[0] -24*param[0]*param[0]*param[0] + 20*param[0]*param[0]-8*param[0]+4);
    curv /= sqrt(9*param[0]*param[0]*param[0]*param[0] -24*param[0]*param[0]*param[0] + 20*param[0]*param[0]-8*param[0]+4);
    curv /= sqrt(9*param[0]*param[0]*param[0]*param[0] -24*param[0]*param[0]*param[0] + 20*param[0]*param[0]-8*param[0]+4);
    return curv;
}

//param[0] = t
void parabola(double*cir, double* param){
    cir[0] = param[0];
    cir[1] = param[0]*param[0];
    cir[2] = 0;
}

void Dparabola(double*p,double *param){
    p[0] = 1;
    p[1] = 2*param[0];
    p[2] = 0;
}

void DDparabola(double*p,double *param){
    if(my_abs(param[0])<1e-5 || my_abs(param[0]-2)<1e-5){
        param[0] +=0.0001;
    }
    p[0] = 0;
    p[1] = 2;
    p[2] = 0;
}

double parabola_curv(double *param){
    // 2/np.sqrt(4*t**2 +1)**3
    double curv = 2;
    curv /= sqrt(4*param[0]*param[0]+1);
    curv /= sqrt(4*param[0]*param[0]+1);
    curv /= sqrt(4*param[0]*param[0]+1);
    return curv;
}