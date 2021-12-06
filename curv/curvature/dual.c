#include<stdio.h>
#include<stdlib.h>
#include "curvature.h" 
#include<math.h>

double my_abs(double x){
    return x>0? x : -x;
}

void printp(double * point){
        printf("[%f %f %f]\n",point[0],point[1],point[2]);
}

double len(double* p1, double* p2){
    return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2] - p2[2])*(p1[2] - p2[2]));
}

void get_point_by_index(double* point,double* points, int indx){
    point[0] = points[3*indx]; point[1] = points[3*indx +1];
    point[2] = points[3*indx +2];
}

void dual_point(double* p1, double* p2, double* n1, double* n2, double* answ){
    double eps = 1E-4;
    //Проверяем нормали на нормированность
    double zeros[3]; zeros[0] = 0; zeros[1] = 0; zeros[2] = 0;
    double l1 = len(n1,zeros);
    double l2 = len(n2,zeros);

    double C1 = -n1[0] * p1[0]/l1 - n1[1]*p1[1]/l1;
    double C2 = -n2[0] * p2[0]/l2 - n2[1]*p2[1]/l2;
    //Случай, когда две нормали коллинеарны
    //Примем тогда дуальную точку - середину
    if(my_abs(n1[1]*n2[0] - n2[1]*n1[0]) < eps ){
        answ[0] = (p1[0] + p2[0])/2;
        answ[1] = (p1[1] + p2[1])/2;
    }else if(my_abs(n1[0])> eps){
        
        answ[1] = ((C1*n2[0]/l2)-(C2*n1[0]/l1))/((n1[0]/l1 )*(n2[1]/l2) - (n2[0]/l2)*(n1[1]/l1));
        answ[0] = (-C1-(n1[1]*answ[1]/l1))/(n1[0]/l1);
    }else if(my_abs(n1[1])> eps){
        answ[0] = ((n2[1]*C1/l2) - (n1[1]*C2/l1))/(((n1[1]*n2[0])/(l1*l2)) - ((n2[1]*n1[0])/(l1*l2)));
        answ[1] = (-C1-(n1[0]*answ[0]/l1))/(n1[1]/l1) ;
    }
    answ[2]=0;
//printf("    answ = ");printp(answ);
}

double* get_dual_points(int* faces, double *points, double* normals,int numFaces){
    //временные переменные, используемые на каждой итерации цикла
    double dualPoint[3],p1[3],p2[3],n1[3],n2[3];
    //Массив, в котором будут хранится старые точки и новые дуальные
    double* dualPoints = (double *) malloc(numFaces*3*sizeof(double));
    //Размер массива с обычными точками
    int i;
    for(i = 0; i<numFaces; i+=1){
        //Из массива точек и  нормалей дастаём текущие точки и нормали
        get_point_by_index(p1,points,faces[3*i]);
        get_point_by_index(p2,points,faces[3*i+1]);
        get_point_by_index(n1,normals,faces[3*i]);
        get_point_by_index(n2,normals, faces[3*i+1]);
        /*printf("i = %d\n",i);
        printf("p1 = ");printp(p1);
        printf("p2 = ");printp(p2);
        printf("n1 = ");printp(n1);
        printf("p2 = ");printp(n2);*/
        //ВЫчисляем дуальную точку по двум точкам и их нормалям
        dual_point(p1,p2, n1, n2, dualPoint);

        //printf("dualPoint = ");printp(dualPoint);
        //Перезаписываем точки
        dualPoints[i*3] = dualPoint[0]; dualPoints[i*3+1] = dualPoint[1];
        dualPoints[i*3+2]=dualPoint[2];    
    }
    //Дополнительная итерация в случае замкнутой кривой
    return dualPoints;
}

int* build_dual_faces(double* points,double* normals,double* dualPoints, int numPoints, int numDualPoints){
    int i,j;
    double point[3],pj[3],p1[3],p2[3],normal[3];
    double curLen, prevLen1, prevLen2,c;
    double dot1,dot2,dot;
    int* dualFaces = (int *) malloc(numPoints* sizeof(int)*3);
    int indxp1 ,indxp2;
    double zero[3];
    double lenn;
    zero[0] = 0 ;zero[1] = 0; zero[2] = 0; 
    for(i = 0;i<numPoints;i++){
        point[0] = points[3*i];point[1] = points[3*i+1]; point[2] = points[3*i+2];
        indxp1 = -1; indxp2 = -1;
        prevLen1 = 1E2; prevLen2 = 2E2; 
        for(j = 0;j<numDualPoints;j++){
            get_point_by_index(pj,dualPoints,j);
            curLen = len(point,pj);
            if(curLen<prevLen2 && curLen<prevLen1){
                prevLen2 = prevLen1;
                prevLen1 = curLen;
                indxp2 = indxp1;
                indxp1 = j;
            }else if(curLen<prevLen2 && curLen>=prevLen1){
                prevLen2 = curLen;
                indxp2 = j;
            }
        }
        if(indxp1 == -1 || indxp2 == -1){
            printf("Something wrong, can't find nearist dual Face for point ");printp(point);
            exit(1);
        }
        get_point_by_index(p1,dualPoints,indxp1);
        get_point_by_index(p2,dualPoints,indxp2);
        /*printf("dualPoint1[%d] = ",indxp1);printp(p1);
        printf("dualPoint2[%d] = ",indxp2);printp(p2);
        printf("for point[%d] ",i); printp(point);*/
        get_point_by_index(normal,normals,i);
        lenn = len(normal,zero); 
        dot1 = (p1[0] - point[0])*normal[0] + (p1[1] - point[1])*normal[1];
        dot2 = (p2[0] - point[0])*normal[0] + (p2[1] - point[1])*normal[1];
        p1[0] = p1[0] - normal[0]*dot1/lenn; p1[1] = p1[1] - normal[1]*dot1/lenn;
        p1[2] = p1[2] - normal[2]*dot1/lenn; p2[0] = p2[0] - normal[0]*dot2/lenn;
        p2[1] = p2[1] - normal[1]*dot2/lenn; p2[2] = p2[2] - normal[2]*dot2/lenn;
        dot = (p1[0] - point[0])*(p2[0] - point[0]) + (p1[1] - point[1])*(p2[1] - point[1]);
        if(dot>0){
            if(numPoints - numDualPoints !=1){
                continue;
            }
            get_point_by_index(p1,dualPoints,indxp1);
            get_point_by_index(p2,dualPoints,indxp2);
            /*printf("p1 = ");printp(p1);
            printf("p2 = ");printp(p2);
            printf("p = "); printp(point);
            printf("normal "); printp(normal);*/
            prevLen2 = 1E2;
            int isFind = 0;
            for(j = 0;j<numDualPoints;j++){
                get_point_by_index(pj,dualPoints,j);
                curLen = len(point,pj);
                dot1 = (p1[0] - point[0])*normal[0] + (p1[1] - point[1])*normal[1];
                dot2 = (pj[0] - point[0])*normal[0] + (pj[1] - point[1])*normal[1];
                p1[0] = p1[0] - normal[0]*dot1/lenn; p1[1] = p1[1] - normal[1]*dot1/lenn;
                p1[2] = p1[2] - normal[2]*dot1/lenn; pj[0] = pj[0] - normal[0]*dot2/lenn;
                pj[1] = pj[1] - normal[1]*dot2/lenn; pj[2] = pj[2] - normal[2]*dot2/lenn;
                dot = (pj[0] - point[0])*(pj[0] - point[0]) + (p1[1] - point[1])*(pj[1] - point[1]);
                if(curLen <prevLen2 && dot<0 ){
                    indxp2 = j;
                    isFind = 1;
                }
            }
            if(isFind == 0){
                 dualFaces[3*i] = -1; dualFaces[3*i + 1] = -1; dualFaces[3*i + 2] = -1;
            }else{
                dualFaces[3*i] = numPoints + indxp1; 
                dualFaces[3*i + 1] =numPoints + indxp2; dualFaces[3*i + 2] =numPoints + indxp2;
            }
        }else{
            dualFaces[3*i] = numPoints + indxp1; 
            dualFaces[3*i + 1] =numPoints + indxp2; dualFaces[3*i + 2] =numPoints + indxp2;
        }
    }
    return dualFaces;
}

int find_dual_face(double* point, double* dualPoints, int* dualFaces, int numDualPoints, int numDualFaces){
    double p1[3],p2[3];
    double curLen, prevLen = 2E3;
    int i;
    int indx = -1;
    for(i = 0; i < numDualFaces; i++){
        if(dualFaces[3*i] == -1) continue;
        get_point_by_index(p1,dualPoints,dualFaces[3*i]);
        get_point_by_index(p2,dualPoints,dualFaces[3*i +1]);
        curLen = len(point,p1) + len(point,p2);
        if(curLen < prevLen){
            indx = i;
            prevLen = curLen;
        }
    }
    return indx;
}

void get_projection(double* projection,double* point, double* p1, double* p2){
    double eps = 1E-7;
    double C,C1,A,B;
    C = p1[0]*(p1[1]-p2[1])+p1[1]*(p2[0]-p1[0]);
    A = p2[1] - p1[1];
    B = p1[0] - p2[0];
    C1 = A*point[1] - B*point[0];
    if(my_abs(p1[0] - point[0])<eps && my_abs(p1[1] - point[1])<eps ){
        projection[0] = point[0];
        projection[1] = point[1];
    } else if(my_abs(p2[0] - point[0])<eps && my_abs(p2[1] - point[1])<eps ){
        projection[0] = point[0];
        projection[1] = point[1];
    }else if(my_abs(B)>eps){
       projection[0] = (-C*A - B*C1)/(A*A+B*B);
       projection[1] = -(C+A*projection[0])/B;
    }else if(my_abs(A)>eps){
       projection[1] = (C1*A - C*B)/(A*A+B*B);
       projection[0] = - (C+B*projection[1])/A;
    }else{
        printf("Something wrong. It seems that points is equal\n");
        printf("A = %f, B = %f\n",A,B);
        printp(p2);
        printp(p1);
        exit(1);
    }
    projection[2] = 0;
}

int find_nearest_face(double* point, double* points, int* faces,int numPoints,int numFaces){
    int j=-1, i;
    double projection[3], p1[3], p2[3];
    double dx,dy,inner, lenn;
    double lenprev = 1E6;
    int istart,iend;
    for(i = 0;i<numFaces;i++){
        if(faces[3*i] == -1) continue;
        get_point_by_index(p1,points,faces[3*i]); 
        get_point_by_index(p2,points,faces[3*i+1]);
        get_projection(projection,point,p1,p2);
        lenn = len(point,projection);
        /*printf("p1 = ");printp(p1);
        printf("p2 = ");printp(p2);
        printf("proj = ");printp(projection);
        printf("len1 = %f, len2 = %f, lenf = %f ,lenn = %f, lenprev = %f \n",len1,len2,lenf,lenn,lenprev);
        printf("len1+len2-lenf = %f \n",my_abs(len1+len2-lenf));*/
        dx = p2[0] - p1[0];
        dy = p2[1] - p1[1];
        inner = (point[0] - p1[0])*dx + (point[1] - p1[1])*dy;
        if( (0 <= inner && inner <= dx*dx + dy*dy) && lenn<lenprev){
            j = i; lenprev = lenn;
        }
    }
    if(j == -1 && numFaces != numPoints){
        for(i = 0; i<numFaces; i++){
            get_point_by_index(p1,points,faces[3*i]); 
            get_point_by_index(p2,points,faces[3*i+1]);
            get_projection(projection,point,p1,p2);
            lenn = len(point,projection);
            if(lenn<lenprev ){
                j = i; lenprev = lenn;
            }
        }
    }
    return j;
}