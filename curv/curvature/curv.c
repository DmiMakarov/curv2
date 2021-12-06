#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "curvature.h"

double get_step(double* p, double* c, double* dc, double* ddc){
    double up = (c[0] - p[0])*dc[0] + (c[1] - p[1])*dc[1] + (c[2] - p[2])*dc[2];
    double d = (p[0] - c[0])* (ddc[0]) + (p[1] - c[1])* (ddc[1]) + (p[2] - c[2])* (ddc[2]); 
    d = d - dc[0]*dc[0] - dc[1]*dc[1] - dc[2]*dc[2];
    return up/d ;
}

/*double get_alpha(double alpha0, double theta, double eps, double step,double curt, void (*fig)(double *,double *), double* dc,double *param){
    double c1[3],c2[3];
    param[0] = curt - alpha0*step;
    fig(c1,param);
    param[0] = curt;
    fig(c2,param);
    while((my_abs(c1[0]- c2[0] + eps*step*alpha0*dc[0]) > eps || my_abs(c1[1] - c2[1] + eps*step*alpha0*dc[1])>eps) || (c1[0] >  c2[0] - eps*step*alpha0*dc[0] || c1[1] > c2[1] < eps*step*alpha0*dc[1] ) && alpha0>=0.0001 ){
        alpha0 = theta*alpha0;
        param[0] = curt - alpha0*step;
        fig(c1,param);
        //printf("    alpha0 = %f, diff1 = %f, diff2 = %f, c[0] = %f, c[1] = %f \n",alpha0, c1[0] - c2[0] + eps*alpha0*dc[0]*step,c1[1]-c2[1] + eps*alpha0*dc[1]*step,c2[0],c1[0]);
        //printf("    theta = %f,eps =%f,step = %f,curt = %f,dc[0] = %f , dc[1]=%f\n",theta,eps,step,curt,dc[0],dc[1]);
    }
    return alpha0;
}*/

double get_progection(double *point, void (*fig)(double *,double *), void (*Dfig)(double*,double*),void (*DDfig)(double *,double*), double* param,double t0,double t1,double tstep,double alpha0,double theta,double eps){
    double c[3],dc[3], ddc[3];
    fig(c,param); Dfig(dc,param); DDfig(ddc,param);
    double a,b;
    double l=1e2;
    int i, num = (int) (t1 - t0)/tstep + 1;
    while((t0 +(num-1)*tstep)<t1) num+=1;
    double curt,tansw;
    double alpha,step;
    int maxnum = 1e5;
    int n;
    for(i = 0;i<num;i++){
        param[0] = t0 + (i-1)*tstep;
        fig(c,param); Dfig(dc,param);
        a = (point[0] - c[0])*dc[0]+(point[1] - c[1])*dc[1]+(point[2] - c[2])*dc[2];
        param[0] = t0 + (i+1)*tstep;
        fig(c,param); Dfig(dc,param);
        b = (point[0] - c[0])*dc[0]+(point[1] - c[1])*dc[1]+(point[2] - c[2])*dc[2];
        //printf("t = %f, a = %f, b = %f\n",t0+i*tstep,a,b);
        if(a*b<=0){
            param[0] = t0 + i*tstep;
            fig(c,param); Dfig(dc,param); DDfig(ddc,param);
            step = get_step(point,c,dc,ddc);
            alpha = alpha0;
            /*if(my_abs(step)>1e-6){
                alpha = get_alpha(alpha0,theta,eps,my_abs(step),param[0],fig,dc,param);
            }*/
            curt = param[0] + alpha*step;
            /*printf("tstart = %f\n",param[0]);
            printf("c = ");printp(c);*/
            n = 0;
            while(my_abs(curt-param[0])>eps && n<maxnum){
            //printf("t1 = %f, t0 = %f\n",t1,param[0]);
                //printf("%f\n",my_abs(curt-param[0]));
                param[0] = curt;
                fig(c,param); Dfig(dc,param); DDfig(ddc,param);
                step = get_step(point,c,dc,ddc);
                /*if(my_abs(step)>eps){
                    alpha = get_alpha(alpha0,theta,eps,my_abs(step),param[0],fig,dc,param);
                }*/
                //printf("alpha = %f\n",alpha);
                curt = param[0] + alpha*step;
                n+=1;
            }
            param[0] = curt;
            fig(c,param);
            //printf("t0 = %f. t1  = %f, ti = %f\n",t0,t1,param[0]);
            if(len(point,c)<l){
                tansw = curt;
                l = len(point,c);
            }
        }
    }
    /*printf("c = ");printp(c);printf("dc = ");printp(dc);printf("ddc = "); printp(ddc);
    printf("t1 = %f, t0 = %f\n",t1,param[0]);*/
    
    return tansw;
}

void get_coef(double* coef,double *p1,double* p2 ){
    coef[0] = p2[1]-p1[1];
    coef[1] = p1[0]-p2[0];
    coef[2] = (p2[0]-p1[0])*p1[1] + (p1[1] - p2[1])*p1[0];
}

double get_angle(double* n1, double* n2){
    double zeros[3]; zeros[0] = 0; zeros[1] = 0; zeros[2] = 0;
    double angle = my_abs(n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2]);
    angle =angle/(len(n1,zeros)*len(n2,zeros)); 
    angle = acos(angle);
    if(angle > pi/2) angle = pi - angle;
    return angle;
}

double get_len(double* a,double *b, double *n1){
    double l = (b[0] - a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1] - a[1]) + (b[2] - a[2])*(b[2] - a[2]);
    double zeros[3]; zeros[0] = 0;zeros[1] = 0; zeros[2] = 0;
    double p = len(n1,zeros);
    double n[3]; 
    n[0] = n1[0]/p; n[1] = n1[1]/p; n[2] = n1[2]/p;
    l =l-(b[0] - a[0])*(b[0]- a[0])*n[0]*n[0];
    l =l-(b[1] - a[1])*(b[1]- a[1])*n[1]*n[1];
    l =l-(b[2] - a[2])*(b[2]- a[2])*n[2]*n[2];
    l =l-2*(b[0] - a[0])*(b[1] - a[1])*n[0]*n[1];
    l =l-2*(b[0] - a[0])*(b[2] - a[2])*n[0]*n[2];
    l =l-2*(b[2] - a[2])*(b[1] - a[1])*n[1]*n[2];
    return sqrt(l);
}

double find_on_curv(double* point,void (*fig)(double *,double *),double t0,double t1,double tstep,double *param){
    double l = 1e2,prevlen = 1e2;
    double c[3];
    int num = (t1-t0)/tstep + 1;
    int i,j;
    while((t0 +(num-1)*tstep)<t1) num+=1;
    for(i = 0;i<num;i++){
        param[0] = t0 + i*tstep;
        fig(c,param);
        l = len(c,point);
        if(l < prevlen){
            j = i;
            prevlen = l;
        }
        if( l < 1e-6) break;
    }
    return t0+j*tstep;
}

void compute_average(double* answ,double *p1,double *p2,double* normal1,double* normal2,double* normalts, double ts,void (*fig)(double *,double *), void (Dfig)(double*,double*),void (DDfig)(double *,double*),double (*figCurv)(double*), double* param,double t0,double t1,double tstep,double alpha0,double theta,double eps){
    double curt1 = find_on_curv(p1,fig,t0,t1,tstep/10,param);
    double curt2 = find_on_curv(p2,fig,t0,t1,tstep/10,param);
    
    int i;
    answ[0] = 0; answ[1] = 0;
    tstep = 0.1*tstep;
    if(curt1 > ts && ts>curt2){
        double c;
        c = curt1;
        curt1 = curt2;
        curt2 = c;
    }
    int num1 = ((ts - curt1)/tstep) + 1;
    int num2 = curt2> curt1 ? ((curt2 - ts)/tstep) + 1 : ((t1-ts)/tstep) + 1;
    double c1[3],c2[3];
    param[0] = num1;
    fig(c1,param);
    param[0] = num2;
    fig(c2,param);
    /*printf("p1 = "); printp(p1);
    printf("p2 = "); printp(p2);
    printf("curt1 = %f, curt2 = %f, ts = %f, num1 = %d, num2 = %d\n",curt1,curt2,ts,num1,num2);
    printf("c1 = ");printp(c1);
    printf("c2 = "); printp(c2);*/
    if(curt1>curt2 && curt1 >ts && ts<curt2){
        //printf("Hello\n");
        num1 = (t1 - curt1)/tstep + ts/tstep +2;
        num2 = ((curt2 - ts)/tstep) + 1;
    }
    for(i= 0;i<num1;i++){
        param[0] = curt1 + i*tstep;
        answ[0] += figCurv(param);
    } 
    for(i = 0;i<num2;i++){
        param[0] = ts + i*tstep;
        answ[1] += figCurv(param);
    }
    answ[0] = answ[0]/num1;
    answ[1] = answ[1]/num2;

    //printf("answ[0] = %f, answ[1] = %f\n",answ[0],answ[1]);
};

curvs* compute_curvs(double *points, double* normals, int *faces, double *dualPoints,int *dualFaces, int numPoints, int numFaces,void (*fig)(double *,double *), void (*Dfig)(double*,double*),void (*DDfig)(double*,double*), double (*figCurv)(double*), double* param,int numFigPoints, double t0, double tmax,int isFlat){
    curvs *curv = (curvs *) malloc(sizeof(curvs));
    curv->w1_ = (double *) malloc(numFaces*sizeof(double));
    curv->w1 = (double *) malloc(numFaces*sizeof(double));
    curv->w2_ = (double *) malloc(numFaces*sizeof(double));
    curv->w2 = (double *) malloc(numFaces*sizeof(double));
    curv->curv = (double *) malloc(numFigPoints*sizeof(double));
    curv->dualCurv = (double *) malloc(numFigPoints*sizeof(double));    
    curv->trueCurv = (double *) malloc(numFigPoints*sizeof(double));
    curv->lens = (double *) malloc(numFigPoints*sizeof(double));
    double *tmpcurv = (double *) malloc(numFaces*2*sizeof(double));
    double *aPoints = (double *) calloc(numFaces*3,sizeof(double)); 
    double *ts = (double *) malloc(numFaces*sizeof(double));
    double tstep = (tmax - t0)/(numFigPoints - 1);
    double phi;
    int i;
    double t;
    int indFace;
    double curPoint[3],c[3];
    double dualPoint[3],point1[3],point2[3];
    double dot,l;
    double dualNormal[3];
    double coef1[3],coef2[3];
    double zeros[3],a[3];
    double normal1[3],normal2[3],acurv[2];
    int isChange = 0;
    zeros[0] = 0; zeros[1] = 0; zeros[2] = 0;
    curv->lens[0] = 0;
    for(i = 0; i<numFigPoints;i++){
        t = t0 + i*tstep;
        param[0] = t;
        fig(curPoint,param);
        if(i!= 0 ){ 
            param[0] = t-tstep;
            fig(c,param);
            curv->lens[i] = curv->lens[i-1] + len(curPoint,c);
            param[0] = t;
        }
        indFace = find_nearest_face(curPoint,points,faces,numPoints,numFaces);
        //printf("indFace = %d \n",indFace);
        if(indFace == -1){
            printf("Can't find suitable faces for point "); printp(curPoint);
        }
        if(aPoints[3*indFace] != 0 || aPoints[3*indFace+1] != 0 ){
            if( i == 0){
                printf("Something wrong in compute_curvs\n");
                break;
            }
            if(t<ts[indFace]){
                curv->curv[i] = curv->w1_[indFace];
                curv->dualCurv[i] = curv->w2_[indFace];
                if(isFlat == 1){
                    curv->trueCurv[i] = tmpcurv[2*indFace];
                }else{
                    param[0] = t;
                    curv->trueCurv[i] = figCurv(param);
                }
            }else{
                curv->curv[i] = curv->w1[indFace];
                curv->dualCurv[i] = curv->w2[indFace];
                if(isFlat == 1){
                    curv->trueCurv[i] = tmpcurv[2*indFace + 1];
                }else{
                    param[0] = t;
                    curv->trueCurv[i] = figCurv(param);
                }
            }
            continue;
        }
        get_point_by_index(dualPoint,dualPoints,indFace);
        get_point_by_index(point1,points,faces[3*indFace]);  
        get_point_by_index(point2,points,faces[3*indFace +1]);
        param[0] = t+tstep;
        fig(c,param);
        param[0] = t;
        dot = (point2[0] - point1[0])*(c[0] - curPoint[0]);
        dot+=(point2[1] - point1[1])*(c[1] - curPoint[1]);
        dot+=(point2[2] - point1[2])*(c[2] - curPoint[2]);
        if(dot<0){
            c[0] = point1[0]; c[1] = point1[1]; c[2] = point1[2];
            point1[0] = point2[0]; point1[1] = point2[1]; point1[2] = point2[2];
            point2[0] = c[0]; point2[1] = c[1]; point2[2] = c[2];
            isChange = 1;
        }
        ts[indFace] = get_progection(dualPoint,fig,Dfig,DDfig,param,t0,tmax,2*tstep,1,0.9,1e-4);
       /* if(ts[indFace]<0 && my_abs(t0)<=1e-3) {
            ts[indFace] = my_abs(ts[indFace]);
        }*/
        if(ts[indFace]>tmax && numFaces == numPoints){
            ts[indFace] = -2*pi + ts[indFace];
        }
        param[0] = ts[indFace];
        fig(c,param);
        get_coef(coef1,dualPoint,c);
        param[0] = t;
        get_coef(coef2,point2,point1);
        aPoints[3*indFace] = (coef2[2]*coef1[1] - coef1[2]*coef2[1])/(coef1[0]*coef2[1] - coef2[0]*coef1[1]);
        aPoints[3*indFace + 1] = (coef1[2]*coef2[0] - coef2[2]*coef1[0])/(coef1[0]*coef2[1] - coef2[0]*coef1[1]);
        a[0] = aPoints[3*indFace];
        a[1] = aPoints[3*indFace+1]; a[2] = 0;
        dualNormal[0] = (dualPoint[0] - a[0]); dualNormal[1] = (dualPoint[1] - a[1]);
        dualNormal[2] = 0;
        l = len(dualNormal,zeros);
        dualNormal[0] = dualNormal[0]/l; dualNormal[1] = dualNormal[1]/l;
        /*if(indFace<5){
        printf("indFace = %d, i = %d, ts = %f\n",indFace,i,ts[indFace]);
        printf("isChange = %d\n",isChange);
        printf("coef1 = [%f,%f,%f]\n",coef1[0],coef1[1],coef1[2]);
        printf("coef2 = [%f,%f,%f]\n",coef2[0],coef2[1],coef2[2]);
        printf("point1 = [%f,%f,%f]\n",point1[0],point1[1],point1[2]);
        printf("point2 = [%f,%f,%f]\n",point2[0],point2[1],point2[2]);
        printf("dualPoint = [%f,%f,%f]\n",dualPoint[0],dualPoint[1],dualPoint[2]);
        printf("c = [%f,%f,%f]\n",c[0],c[1],c[2]);
        printf("a = "); printp(a);
        printf("lenDual = %f\n",len(dualNormal,zeros));
        }*/
        if(isChange == 0){
            get_point_by_index(normal1,normals,faces[3*indFace]);
            get_point_by_index(normal2,normals,faces[3*indFace + 1]);
        }else{
            get_point_by_index(normal1,normals,faces[3*indFace+1]);
            get_point_by_index(normal2,normals,faces[3*indFace]);
        }
        phi = get_angle(normal1,dualNormal);
        /*if(indFace < 5){
        printf("for phi1 = %f , get_len = %f, len = %f,normal1 ",phi,get_len(a,point1,dualNormal),get_len(point1,dualPoint,normal1));printp(normal1);
        }*/
        curv->w1_[indFace] = phi/get_len(a,point1,dualNormal);
        curv->w2_[indFace] = phi/get_len(point1,dualPoint,normal1);
        phi = get_angle(normal2,dualNormal);
        /*if(indFace<5){
        printf("for phi2 = %f , get_len = %f, len = %f,normal2 ",phi,get_len(a,point2,dualNormal),get_len(point2,dualPoint,normal2));printp(normal2);
        printf("dualNormal "); printp(dualNormal);
        }*/
        curv->w1[indFace] = phi/get_len(a,point2,dualNormal);
        curv->w2[indFace] = phi/get_len(point2,dualPoint,normal2);
        /*if(indFace<5){
        printf("w1_[%d] = %f, w1[%d] = %f\n",indFace,curv->w1_[indFace],indFace,curv->w1[indFace]);
        printf("w2_[%d] = %f, w2[%d] = %f\n",indFace,curv->w2_[indFace],indFace,curv->w2[indFace]);
        }*/
        if(isFlat == 1){
            compute_average(acurv,point1,point2,normal1,normal2,dualNormal,ts[indFace],fig,Dfig,DDfig,figCurv,param,t0,tmax,tstep,1,0.9,1e-5);
            tmpcurv[2*indFace] = acurv[0]; tmpcurv[2*indFace+1] = acurv[1];
            if(t<ts[indFace]){
                curv->curv[i] = curv->w1_[indFace];
                curv->dualCurv[i] = curv->w2_[indFace];
                curv->trueCurv[i] = acurv[0];
            }else{
                curv->curv[i] = curv->w1[indFace];
                curv->dualCurv[i] = curv->w2[indFace];
                curv->trueCurv[i] = acurv[1];
            }
        /*printf("param = [%f,%f]\n",param[0],param[1]);
        printf("curv = %f\n",figCurv(param));*/
        }else{
            param[0] = t;
            curv->trueCurv[i] = figCurv(param); 
            if(t<ts[indFace]){
                curv->curv[i] = curv->w1_[indFace];
                curv->dualCurv[i] = curv->w2_[indFace];
            }else{
                curv->curv[i] = curv->w1[indFace];
                curv->dualCurv[i] = curv->w2[indFace];
            }
        }
        //printf("curv = %f, dualCurv = %f, truecurv = %f\n",curv->curv[i],curv->dualCurv[i],curv->trueCurv[i]);
        isChange = 0;
    }
    free(aPoints);
    free(ts);
    return curv;
}

void free_curvs(curvs* curv){
    free(curv->w1); free(curv->w1_); free(curv->w2);
    free(curv->w2_); free(curv->curv); free(curv->dualCurv);
    free(curv->lens); free(curv->trueCurv); free(curv);
}