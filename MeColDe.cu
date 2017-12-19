// Reference Scan implementation - Author: Ananoymous student of ME759 Fall 2017
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include <cuda.h>
#include <math.h>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/scan.h>


void initializeSize(FILE* fp, int* nVertices, int* nFaces)
{
	char x[1024];
	while (fscanf(fp, "%1023s", x)) {
        	if (strcmp(x,"ply")==0) continue;
		if ((strcmp(x,"format")==0)||(strcmp(x,"comment")==0)||(strcmp(x,"property")==0)) {fgets(x, 1023, fp);continue;}
								//	{fscanf(fp,"%[^\n]", x);continue;}
		if (strcmp(x,"element")==0) {
			fscanf(fp, "%1023s", x);
			if (strcmp(x,"vertex")==0) fscanf(fp, "%d\n", &nVertices[0]);
			else if (strcmp(x,"face")==0) fscanf(fp, "%d\n", &nFaces[0]);
		}
		if (strcmp(x,"end_header")==0) break;
    	}
}

void initializeArr(FILE* fp, double* ArrVerX, double* ArrVerY, double* ArrVerZ, int* ArrFace, int V, int F)
{
	int i;
	char x[1024];	
	for(i=0; i<V; ++i){
		fscanf(fp,"%f",&ArrVerX[i]);
		fscanf(fp,"%f",&ArrVerY[i]);
		fscanf(fp,"%f",&ArrVerZ[i]);
		fgets(x, 1023, fp);
//		if(r == EOF){
//			rewind(fp);
//		}
	}

	for (i=0; i<F; ++i){
		fscanf(fp, "%1023s", x);
		fscanf(fp, "%d", &ArrFace[3*i]);
		fscanf(fp, "%d", &ArrFace[3*i+1]);
		fscanf(fp, "%d", &ArrFace[3*i+2]);
	}
}

bool intersect(int i, int j, double *p_x, double *p_y, double *p_z, double *tri)
{



}

double get_volume(int num_1, int num_2, double* faces_1, double* faces_2)
{


}

int main(int argc, char* argv[]) {
	FILE *fp_1 = fopen("bun_zipper.ply","r");
	FILE *fp_2 = fopen("test_mesh.ply","r");
	//allocate resources
	int nV_1=0, nF_1=0, nV_2=0, nF_2=0,i,j;
	int num_ray_intersect;
	float time = 0.f;
	initializeSize(fp_1,&nV_1, &nF_1);
	initializeSize(fp_2,&nV_2, &nF_2);
	double *vertices_x_1= (double *)malloc(sizeof(double)*nV_1);
	double *vertices_y_1= (double *)malloc(sizeof(double)*nV_1);
	double *vertices_z_1= (double *)malloc(sizeof(double)*nV_1);
        int *faces_1   = (int *)malloc(sizeof(int)*nF_1*3);
	double *vertices_x_2= (double *)malloc(sizeof(double)*nV_2);
        double *vertices_y_2= (double *)malloc(sizeof(double)*nV_2);
        double *vertices_z_2= (double *)malloc(sizeof(double)*nV_2);
        int *faces_2   = (int *)malloc(sizeof(int)*nF_2*3);

	initializeArr(fp_1,vertices_x_1,vertices_y_1,vertices_z_1,faces_1,nV_1,nF_1);
	initializeArr(fp_2,vertices_x_2,vertices_y_2,vertices_z_2,faces_2,nV_2,nF_2);
	// Your code here
 
  	double *face_coord_1 = (double *)malloc(sizeof(double)*nF_1*9);
	double *face_coord_2 = (double *)malloc(sizeof(double)*nF_2*9);
 // cudaMalloc((void**)&dout,size);
	//cudaMalloc((void**)&din,size);
	int *inside_point_set_1 = (int *)malloc(sizeof(int)*nV_1);
	int *inside_point_set_2 = (int *)malloc(sizeof(int)*nV_2);
 
  	cudaEvent_t startEvent_inc, stopEvent_inc;
	cudaEventCreate(&startEvent_inc);
	cudaEventCreate(&stopEvent_inc);
  	cudaEventRecord(startEvent_inc,0); // starting timing for inclusive  

	for (i=0;i<nF_1;++i){
		face_coord_1[9*i] = vertices_x_1[faces_1[3*i]];
		face_coord_1[9*i+1] = vertices_y_1[faces_1[3*i]];
  		face_coord_1[9*i+2] = vertices_z_1[faces_1[3*i]];
		face_coord_1[9*i+3] = vertices_x_1[faces_1[3*i+1]];
                face_coord_1[9*i+4] = vertices_y_1[faces_1[3*i+1]];
                face_coord_1[9*i+5] = vertices_z_1[faces_1[3*i+1]];
		face_coord_1[9*i+6] = vertices_x_1[faces_1[3*i+2]];
                face_coord_1[9*i+7] = vertices_y_1[faces_1[3*i+2]];
                face_coord_1[9*i+8] = vertices_z_1[faces_1[3*i+2]];
	}

	for (i=0;i<nF_2;++i){
                face_coord_2[9*i] = vertices_x_2[faces_2[3*i]];
                face_coord_2[9*i+1] = vertices_y_2[faces_2[3*i]];
                face_coord_2[9*i+2] = vertices_z_2[faces_2[3*i]];
                face_coord_2[9*i+3] = vertices_x_2[faces_2[3*i+1]];
                face_coord_2[9*i+4] = vertices_y_2[faces_2[3*i+1]];
                face_coord_2[9*i+5] = vertices_z_2[faces_2[3*i+1]];
                face_coord_2[9*i+6] = vertices_x_2[faces_2[3*i+2]];
                face_coord_2[9*i+7] = vertices_y_2[faces_2[3*i+2]];
                face_coord_2[9*i+8] = vertices_z_2[faces_2[3*i+2]];
        }
	//AABB construction
	//NO, dont do this; instead, find points if inside the mesh to determine if this point belongs to the intersection part, this is not likely to be wrong in such a convex struct
	//how to determine if a point is inside the mesh? a lite way, find the closet point on the mesh(para), then move the point towards the direction defined by the point the center of the mesh(the dist to move is around 1/10 of the element length), if the closet point becomes even closer, it should be outside; or it should be inside
	//you have all the insided tris, then ...
	//divergence theorem to locate the volume

	int num_point_inside_1 = 0, num_point_inside_2 = 0;  
	int num_relevant_face_1 = 0, num_relevant_face_2 = 0;
	for (i = 0;i<nV_2;++i) {   //index of points that need to know if inside mesh
		num_ray_intersect = 0;
		for (j=0;j<nV_1;++j) {  //index of triangles
			if (((vertices_y_2[i]>face_coord_1[9*j+1])&&(vertices_y_2[i]>face_coord_1[9*j+4])&&(vertices_y_2[i]>face_coord_1[9*j+7])) ||
			   ((vertices_y_2[i]<face_coord_1[9*j+1])&&(vertices_y_2[i]<face_coord_1[9*j+4])&&(vertices_y_2[i]<face_coord_1[9*j+7])) ||
			   ((vertices_z_2[i]>face_coord_1[9*j+2])&&(vertices_y_2[i]>face_coord_1[9*j+5])&&(vertices_y_2[i]>face_coord_1[9*j+8])) ||
                           ((vertices_z_2[i]<face_coord_1[9*j+2])&&(vertices_y_2[i]<face_coord_1[9*j+5])&&(vertices_y_2[i]<face_coord_1[9*j+8]))) {
				continue;
			}
			else if (intersect(i, j, vertices_x_2, vertices_y_2, vertices_z_2, face_coord_1)) {
				num_ray_intersect++;
			}
		}
		if (num_ray_intersect%2==1) { inside_point_set_2[num_point_inside_2] = i; num_point_inside_2++; }
	}

	for (i = 0;i<nV_1;++i) {   //index of points that need to know if inside mesh
                num_ray_intersect = 0;
                for (j=0;j<nV_2;++j) {  //index of triangles
                        if (((vertices_y_1[i]>face_coord_2[9*j+1])&&(vertices_y_1[i]>face_coord_2[9*j+4])&&(vertices_y_1[i]>face_coord_2[9*j+7])) ||
                           ((vertices_y_1[i]<face_coord_2[9*j+1])&&(vertices_y_1[i]<face_coord_2[9*j+4])&&(vertices_y_1[i]<face_coord_2[9*j+7])) ||
                           ((vertices_z_1[i]>face_coord_2[9*j+2])&&(vertices_y_1[i]>face_coord_2[9*j+5])&&(vertices_y_1[i]>face_coord_2[9*j+8])) ||
                           ((vertices_z_1[i]<face_coord_2[9*j+2])&&(vertices_y_1[i]<face_coord_2[9*j+5])&&(vertices_y_1[i]<face_coord_2[9*j+8]))) {
                                continue;
                        }
                        else if (intersect(i, j, vertices_x_1, vertices_y_1, vertices_z_1, face_coord_2)) {
                                num_ray_intersect++;
                        }
                }
                if (num_ray_intersect%2==1) { inside_point_set_1[num_point_inside_1] = i; num_point_inside_1++; }
        }

	double *relevant_face_1 = (double *)malloc(sizeof(double)*nF_1*9);
        double *relevant_face_2 = (double *)malloc(sizeof(double)*nF_2*9);

	for (i=0;i<nF_1;++i) {
		for (j=0;j<num_point_inside_1;++j) {	
			if ((inside_point_set_1[j]==faces_1[3*i])||(inside_point_set_1[j]==faces_1[3*i+1])||(inside_point_set_1[j]==faces_1[3*i+2])) {
				memcpy(&(relevant_face_1[9*num_relevant_face_1]), &(face_coord_1[9*i]), sizeof(double)*9);
				num_relevant_face_1++;
				break;
			}
		}
	}

	
	for (i=0;i<nF_2;++i) {
                for (j=0;j<num_point_inside_2;++j) {
                        if ((inside_point_set_2[j]==faces_2[3*i])||(inside_point_set_2[j]==faces_2[3*i+1])||(inside_point_set_2[j]==faces_2[3*i+2])) {
                                memcpy(&(relevant_face_2[9*num_relevant_face_2]), &(face_coord_2[9*i]), sizeof(double)*9);
                                num_relevant_face_2++;
                                break;
                        }
                }
        }

	double vol = get_volume(num_relevant_face_1, num_relevant_face_2, relevant_face_1, relevant_face_2);

  	cudaEventRecord(stopEvent_inc,0);  //ending timing for inclusive
  	cudaEventSynchronize(stopEvent_inc);   
	cudaEventElapsedTime(&time, startEvent_inc, stopEvent_inc);   
 

	printf("%d\n%d\n", nV_1, nF_1);


	//free resources 
//	free(in); free(out); free(cuda_out);
	return 0;
}
