#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "curvature.h"
#include "../rw_mesh/rw_mesh.h"

void add_to_grid(struct RW_MESH_VTK_STRUCT *VTK,struct RW_MESH_VTK_UNSTRUCTURED_GRID_STRUCT * Grid, double *dualPoints, int *dualFaces){
    int i;
    int countOfCells = Grid->CountOfCells;
    int countOfPoints = Grid->CountOfPoints;
    int isClosed = countOfCells == countOfPoints ? 1 : 0;
    int sizeNewPoints = isClosed == 1 ? 2*countOfPoints : 2*countOfPoints - 1;
    int sizeNewCells = isClosed == 1 ? 2*countOfCells : 2*countOfCells - 1;
    double* newPoints = (double *) malloc(sizeNewPoints*3*sizeof(double));
    int* newFaces = (int *) malloc(sizeNewCells*3*sizeof(int));
    int* newoffset = (int *) malloc((sizeNewCells+1)*sizeof(int));
    int* newcellsizes = (int*) malloc(sizeNewCells*sizeof(int));
    int* newtypes = (int*) malloc(sizeNewCells*sizeof(int));
    for(i = 0; i<countOfCells; i++){
        newFaces[3*i] = Grid->Cells[3*i];
        newFaces[3*i+1] = Grid->Cells[3*i+1];
        newFaces[3*i+2] = Grid->Cells[3*i+2]; 
    }
    for(i = 0; i<countOfPoints;i++){
        newPoints[3*i] = Grid->Points[3*i];
        newPoints[3*i+1] = Grid->Points[3*i+1];
        newPoints[3*i+2] = Grid->Points[3*i+2];
    }
    for(i = countOfPoints; i<sizeNewPoints; i++){
        newPoints[3*i] = dualPoints[3*(i - countOfPoints)];
        newPoints[3*i+1] = dualPoints[3*(i - countOfPoints)+1];
        newPoints[3*i+2] = dualPoints[3*(i - countOfPoints)+2];
    }
    int j=0;
    for(i = countOfCells; i<sizeNewCells;i++){
        if(dualFaces[3*(j)] == -1){
            j +=1;
            i -=1;
            continue;
        }
        //printf("i = %d,j = %d\n",i,j);
        //printf("[%d,%d,%d]\n",dualFaces[3*j],dualFaces[3*j+1],dualFaces[3*j+2]);
        newFaces[3*i] = dualFaces[3*(j)];
        newFaces[3*i+1] = dualFaces[3*(j)+1];
        newFaces[3*i+2] = dualFaces[3*(j)+2];
        j+=1;
    }
    newoffset[0] = 0;
    for(i = 0; i<sizeNewCells; i++) newcellsizes[i] = 3;
    free(Grid->CellSizes);
    Grid->CellSizes = newcellsizes;
    for(i = 0; i<sizeNewCells; i++){ newoffset[i+1] =newoffset[i]+Grid->CellSizes[i];}
    for(i = 0; i<sizeNewCells; i++){ newtypes[i] = 5;}
    free(Grid->CellOffset);
    free(Grid->CellTypes);
    Grid->CellOffset = newoffset;
    Grid->CellTypes = newtypes;
    free(Grid->Points);
    Grid->Points = newPoints;
    Grid->CountOfPoints = sizeNewPoints;
    free(Grid->Cells);
    Grid->Cells = newFaces;
    Grid->CountOfCells = sizeNewCells;
    
    double *newnormals = (double*) calloc(Grid->CountOfPoints*3,sizeof(double));
    double *newcurv = (double *) calloc(Grid->CountOfPoints,sizeof(double));
    struct RW_MESH_VTK_DATASET_STRUCT * PointsData = (struct RW_MESH_VTK_DATASET_STRUCT* )VTK->PointsData;
    for(i = 0;i<VTK->CountOfPointsData; i++){
        if(strcmp((PointsData+i)->DataName ,"normals") == 0){
            struct RW_MESH_VTK_DATA_NORMALS_STRUCT* normals_data = (struct RW_MESH_VTK_DATA_NORMALS_STRUCT *) (PointsData+i)->Data;
            for(j = 0; j<3*countOfPoints;j++){
                newnormals[j] = normals_data->values[j]; 
            } 
            free(normals_data->values);
            normals_data->values = newnormals;
            (PointsData+i)->Counts = sizeNewPoints;
        }
        if(strcmp((PointsData+i)->DataName ,"curvature") == 0){
            struct RW_MESH_VTK_DATA_NORMALS_STRUCT* curv_data = (struct RW_MESH_VTK_DATA_NORMALS_STRUCT *) (PointsData+i)->Data;
            for(j = 0; j<countOfPoints;j++){
                newcurv[j] = curv_data->values[j]; 
            } 
            free(curv_data->values);
            curv_data->values = newcurv;   
            (PointsData+i)->Counts = sizeNewPoints; 
        }
    }
}

int write_file(curvs *curv,int numPoints, char* filename){
    FILE *fd = fopen(filename,"w");
    int i;
    if (fd == NULL){
        printf("Can't open file %s\n",filename);
        fclose(fd);
        return 1;
    }
    fprintf(fd,"Len trueCurv Curv DualCurv2\n");
    for(i = 0; i<numPoints;i++){
        fprintf(fd,"%f %f %f %f\n",curv->lens[i],curv->trueCurv[i], curv->curv[i], curv->dualCurv[i]);
    }
    fclose(fd);
    return 0;
}

int main(int argc, char* argv[]){
    if(argc <2){
        printf("Invalid argument\n");
        return 0;
    }
    int i,iexist = 0, oexist = 0, interpol = 0, wcurv = 0;
    char ifilename[128],ofilename[128],ocurv[128],figName[16],path[128];
    for(i = 1; i<argc;i++){
        if(strcmp(argv[i], "-i") == 0){
            iexist = 1;
            strcpy(ifilename,argv[i+1]);
        }
        if(strcmp(argv[i],"-o") == 0){
            oexist = 1;
            strcpy(ofilename,argv[i+1]);
        }
        if(strcmp(argv[i],"-f") == 0){
            strcpy(figName,argv[i+1]);
            strcpy(path,argv[i+2]);
        }
    }
    if(iexist != 1 || oexist != 1){
        printf("Invalid argument\n");
        return 0;
    }
    //printf("%s\n",figName);
    struct RW_MESH_VTK_STRUCT VTK;
    struct RW_MESH_VTK_UNSTRUCTURED_GRID_STRUCT *Grid;
    int flagRead = read_format_vtk_struct(&VTK, ifilename, 0);
    if(flagRead != 0 ){
        printf("Can't read input file\n");
        rw_mesh_print_error();
        rw_mesh_vtk_struct_free(&VTK); 
        return 0;
    }
    Grid = (struct RW_MESH_VTK_UNSTRUCTURED_GRID_STRUCT *)VTK.Grid;
    int numPoints = Grid->CountOfPoints;
    int numFaces = Grid->CountOfCells;
    struct RW_MESH_VTK_DATASET_STRUCT * PointsData = (struct RW_MESH_VTK_DATASET_STRUCT* )VTK.PointsData;
    struct RW_MESH_VTK_DATA_NORMALS_STRUCT* normals_data;
    for(i = 0;i<VTK.CountOfPointsData; i++){
        if(strcmp((PointsData+i)->DataName,"normals") == 0){
            normals_data = (struct RW_MESH_VTK_DATA_NORMALS_STRUCT *) (PointsData+i)->Data;
        }
    }
    int isClosed = numFaces == numPoints ? 1 : 0;
    double *dualPoints = get_dual_points(Grid->Cells,Grid->Points,normals_data->values,numFaces);
    
    int* dualFaces = build_dual_faces(Grid->Points,normals_data->values,dualPoints,numPoints,numFaces);
    /*for(i=0;i<numPoints;i++){
        printf("[%d,%d,%d]\n",dualFaces[3*i],dualFaces[3*i+1],dualFaces[3*i+2]);
    }*/
    curvs *curv;
    int isFlat = 1;
    //printf("Stop computing duals\n");
    if(strcmp(figName,"circle") == 0){
        double param[2];
        param[0] = 0; param[1] = 10;
        curv = compute_curvs(Grid->Points, normals_data->values,Grid->Cells,dualPoints, dualFaces,numPoints,numFaces,circle,Dcircle,DDcircle,circle_curv,param,250,0,2*pi, isFlat); 
    }else if(strcmp(figName,"ellipse") == 0){
        double param[3];
        param[0] = 0; param[1] = 10; param[2] = 4;
        curv = compute_curvs(Grid->Points, normals_data->values,Grid->Cells,dualPoints, dualFaces,numPoints,numFaces,ellipse,Dellipse,DDellipse,ellipse_curv,param,250,0,2*pi - 0.001, isFlat);
    }else if(strcmp(figName,"spiral") == 0){
        double param[2];
        param[0] = 0; param[1] = 10;
        curv = compute_curvs(Grid->Points, normals_data->values,Grid->Cells,dualPoints, dualFaces,numPoints,numFaces,spiral,Dspiral,DDspiral,spiral_curv,param,250,0.0001,pi, isFlat);
    }else if(strcmp(figName,"parabola") == 0){
        double param[1];
        param[0] = -5;
        curv = compute_curvs(Grid->Points, normals_data->values,Grid->Cells,dualPoints, dualFaces,numPoints,numFaces,parabola,Dparabola,DDparabola,parabola_curv,param,250,-5,5, isFlat);
    }else if(strcmp(figName,"loop") == 0){
        double param[1];
        param[0] = -1; 
        curv = compute_curvs(Grid->Points, normals_data->values,Grid->Cells,dualPoints, dualFaces,numPoints,numFaces,loop,Dloop,DDloop,loop_curv,param,250,-1,3, isFlat);
    }else if(strcmp(figName,"delta") == 0){
        double param[3];
        param[0] = 0; param[1] = 4; param[2] = 12;
        printf("Hello from delta\n");
        curv = compute_curvs(Grid->Points, normals_data->values,Grid->Cells,dualPoints, dualFaces,numPoints,numFaces,delta,Ddelta,DDdelta,delta_curv,param,250,0.0001,2*pi, isFlat);
    }else{
        printf("Unrnown figure name\n");
        rw_mesh_vtk_struct_free(&VTK);
        free(dualPoints);
        free(dualFaces);
        return 0;
    }
   //printf("Stop computing\n");
    add_to_grid(&VTK, Grid,dualPoints,dualFaces);
    int ret =write_format_vtk_struct(&VTK, ofilename, 0);
    if(ret != 0){
        printf("Can't write VTK\n");
    }
    ret = write_file(curv,250,path);
    if(ret!=0){
        printf("Can't write file\n");
    }
    rw_mesh_vtk_struct_free(&VTK);
    free_curvs(curv);
    free(dualPoints);
    free(dualFaces);
}