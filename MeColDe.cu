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
			else {fgets(x, 1023, fp);continue;}
		}
		if (strcmp(x,"end_header")==0) break;
    	}
}

void initializeArr(FILE* fp, double* ArrVerX, double* ArrVerY, double* ArrVerZ, int* ArrFace, int V, int F)
{
	int i;
	char x[1024];	
	for(i=0; i<V; ++i){
		fscanf(fp,"%lf",&ArrVerX[i]);
		fscanf(fp,"%lf",&ArrVerY[i]);
		fscanf(fp,"%lf",&ArrVerZ[i]);
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

void cross(double* res, double* u, double* v) 
{
	res[0] = u[1]*v[2] - u[2]*v[1];
	res[1] = u[2]*v[0] - u[0]*v[2];
	res[2] = u[0]*v[1] - u[1]*v[0];
}

/*
__device__ void d_cross(double* res, double* u, double* v)
{
        res[0] = u[1]*v[2] - u[2]*v[1];
        res[1] = u[2]*v[0] - u[0]*v[2];
        res[2] = u[0]*v[1] - u[1]*v[0];
}
*/



double dot(double* u, double* v)
{
	return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}

bool intersect(double *p_x, double *p_y, double *p_z, double *tri)
{
	double small_num = 0.00000001;
	double e1[3] = {tri[3]-tri[0], tri[4]-tri[1], tri[5]-tri[2]};
	double e2[3] = {tri[6]-tri[0], tri[7]-tri[1], tri[8]-tri[2]};

//	double norm = sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
//	if (norm == 0.0) {printf("oooh");}

	double d[3] = {0,1,0};
	double h[3];
	cross(h, d, e2);
	double a = dot(e1, h);
	if ((a>(-1*small_num))&&(a<small_num)) return false;
	double f = 1/a;
	double s[3] = {*p_x-tri[0], *p_y-tri[1], *p_z-tri[2]};
	double u = f*(dot(s, h));
	if ((u<0.0)||(u>1.0)) return false;
	double q[3];
	cross(q, s, e1);
	double v = f*(dot(d, q));
	if ((v<0.0)||(u+v>1.0)) return false;
	double t = f*(dot(e2, q));
	if (t>small_num) return true;
	else return false;
}

__global__ void find_volume_normal(int tot_num, double *faces, double *volume, double *normal) {
    	int yourID = blockIdx.x*blockDim.x+threadIdx.x;
	double *yourTri;
	yourTri = faces + 9*yourID;
	if (yourID<tot_num) {
		double d13[3] = {yourTri[3]-yourTri[6], yourTri[4]-yourTri[7], yourTri[5]-yourTri[8]};
		double d12[3] = {yourTri[0]-yourTri[3], yourTri[1]-yourTri[4], yourTri[2]-yourTri[5]};
		double cr[3];
		//d_cross(cr, d13, d12);
		cr[0] = d13[1]*d12[2] - d13[2]*d12[1];
        	cr[1] = d13[2]*d12[0] - d13[0]*d12[2];
        	cr[2] = d13[0]*d12[1] - d13[1]*d12[0];
		
		double crNorm = sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);

//		if (crNorm==0.0) crNorm = 1;	

		double area = 0.5*crNorm;
		double zMean = (yourTri[2]+yourTri[5]+yourTri[8])/3;
		double nz = (-1)*cr[2]/crNorm;
		volume[yourID] = area*zMean*nz;
		normal[3*yourID] = cr[0]/crNorm;
		normal[3*yourID+1] = cr[1]/crNorm;
		normal[3*yourID+2] = cr[2]/crNorm;
	}
}

double get_volume(int num_1, int num_2, double* faces_1, double* faces_2, double* normal_1)
{

	double *d_faces;  
	cudaMalloc(&d_faces, sizeof(double)*(num_1+num_2)*9);

	double *d_vol, *d_nor;
	cudaMalloc(&d_vol, sizeof(double)*(num_1+num_2)*1);
	cudaMalloc(&d_nor, sizeof(double)*(num_1+num_2)*3);

	double *h_vol = (double *)malloc(sizeof(double)*(num_1+num_2)*1);
	double *h_nor = (double *)malloc(sizeof(double)*(num_1+num_2)*3);

	cudaMemcpy(d_faces,faces_1,sizeof(double)*num_1*9,cudaMemcpyHostToDevice);
	cudaMemcpy(d_faces+num_1*9,faces_2,sizeof(double)*num_2*9,cudaMemcpyHostToDevice);
  	find_volume_normal<<<(num_1+num_2+1023)/1024, 1024>>>(num_1+num_2, d_faces, d_vol, d_nor);
  	cudaMemcpy(h_vol, d_vol, sizeof(double)*(num_1+num_2), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_nor, d_nor, sizeof(double)*(num_1+num_2)*3,  cudaMemcpyDeviceToHost);

	double sum_v = 0.0, nor_x_tot = 0.0, nor_y_tot = 0.0, nor_z_tot = 0.0;
	int i;
	for (i = 0; i<num_1+num_2; ++i) {
		sum_v = sum_v + h_vol[i];
	}

	for (i = 0; i<num_1; ++i) {
		nor_x_tot += h_nor[3*i];
		nor_y_tot += h_nor[3*i+1];
		nor_z_tot += h_nor[3*i+2];
	}

	normal_1[0] = nor_x_tot/num_1;
	normal_1[1] = nor_y_tot/num_1;
	normal_1[2] = nor_z_tot/num_1;
	return sum_v;
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

// for test reasons
//	double dis = 0.05;
//	for (i=0;i<nV_2;++i) {
//		vertices_x_2[i] += dis;
//		vertices_y_2[i] += dis;
//		vertices_z_2[i] += dis;
//	}

//	printf("%f\n%d\n", *vertices_x_1, *faces_2);

//end of test

 
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
		for (j=0;j<nF_1;++j) {  //index of triangles
			if (((vertices_x_2[i]>face_coord_1[9*j+0])&&(vertices_x_2[i]>face_coord_1[9*j+3])&&(vertices_x_2[i]>face_coord_1[9*j+6])) ||
			   ((vertices_x_2[i]<face_coord_1[9*j+0])&&(vertices_x_2[i]<face_coord_1[9*j+3])&&(vertices_x_2[i]<face_coord_1[9*j+6])) ||
			   ((vertices_z_2[i]>face_coord_1[9*j+2])&&(vertices_z_2[i]>face_coord_1[9*j+5])&&(vertices_z_2[i]>face_coord_1[9*j+8])) ||
                           ((vertices_z_2[i]<face_coord_1[9*j+2])&&(vertices_z_2[i]<face_coord_1[9*j+5])&&(vertices_z_2[i]<face_coord_1[9*j+8]))) {
				continue;
			}
			else if (intersect(vertices_x_2+i, vertices_y_2+i, vertices_z_2+i, face_coord_1+9*j)) {
				num_ray_intersect++;
			}
		}
		if (num_ray_intersect%2==1) { inside_point_set_2[num_point_inside_2] = i; num_point_inside_2++; }
	}

	for (i = 0;i<nV_1;++i) {   //index of points that need to know if inside mesh
                num_ray_intersect = 0;
                for (j=0;j<nF_2;++j) {  //index of triangles
                        if (((vertices_x_1[i]>face_coord_2[9*j+0])&&(vertices_x_1[i]>face_coord_2[9*j+3])&&(vertices_x_1[i]>face_coord_2[9*j+6])) ||
                           ((vertices_x_1[i]<face_coord_2[9*j+0])&&(vertices_x_1[i]<face_coord_2[9*j+3])&&(vertices_x_1[i]<face_coord_2[9*j+6])) ||
                           ((vertices_z_1[i]>face_coord_2[9*j+2])&&(vertices_z_1[i]>face_coord_2[9*j+5])&&(vertices_z_1[i]>face_coord_2[9*j+8])) ||
                           ((vertices_z_1[i]<face_coord_2[9*j+2])&&(vertices_z_1[i]<face_coord_2[9*j+5])&&(vertices_z_1[i]<face_coord_2[9*j+8]))) {
                                continue;
                        }
                        else if (intersect(vertices_x_1+i, vertices_y_1+i, vertices_z_1+i, face_coord_2+9*j)) {
                                num_ray_intersect++;
                        }
                }
                if (num_ray_intersect%2==1) { inside_point_set_1[num_point_inside_1] = i; num_point_inside_1++; }
        }


	double x_tot = 0.0, y_tot = 0.0, z_tot = 0.0;
	for (i=0;i<num_point_inside_1;++i) {
		x_tot += vertices_x_1[inside_point_set_1[i]];
		y_tot += vertices_y_1[inside_point_set_1[i]];
		z_tot += vertices_y_1[inside_point_set_1[i]];
	}

	for (i=0;i<num_point_inside_2;++i) {
                x_tot += vertices_x_2[inside_point_set_2[i]];
                y_tot += vertices_y_2[inside_point_set_2[i]];
                z_tot += vertices_z_2[inside_point_set_2[i]];
        }

	double center[3] = {x_tot/(num_point_inside_1+num_point_inside_2), 
				y_tot/(num_point_inside_1+num_point_inside_2), 
				z_tot/(num_point_inside_1+num_point_inside_2)};

	double *relevant_face_1 = (double *)malloc(sizeof(double)*nF_1*9);
        double *relevant_face_2 = (double *)malloc(sizeof(double)*nF_2*9);

	for (i=0;i<nF_1;++i) {
		for (j=0;j<num_point_inside_1;++j) {	
			if ((inside_point_set_1[j]==faces_1[3*i])||(inside_point_set_1[j]==faces_1[3*i+1])||(inside_point_set_1[j]==faces_1[3*i+2])) {
				memcpy(relevant_face_1+9*num_relevant_face_1, face_coord_1+9*i, sizeof(double)*9);
				num_relevant_face_1++;
				break;
			}
		}
	}

	
	for (i=0;i<nF_2;++i) {
                for (j=0;j<num_point_inside_2;++j) {
                        if ((inside_point_set_2[j]==faces_2[3*i])||(inside_point_set_2[j]==faces_2[3*i+1])||(inside_point_set_2[j]==faces_2[3*i+2])) {
                                memcpy(relevant_face_2+9*num_relevant_face_2, face_coord_2+9*i, sizeof(double)*9);
                                num_relevant_face_2++;
                                break;
                        }
                }
        }
	
	printf("%d   %d\n",num_point_inside_1, num_point_inside_2);
	printf("%d   %d\n",num_relevant_face_1, num_relevant_face_2);
	printf("%f\n", *relevant_face_1); 	


	double *normal_1 = (double *)malloc(sizeof(double)*3);

	double vol = get_volume(num_relevant_face_1, num_relevant_face_2, relevant_face_1, relevant_face_2, normal_1);

  	cudaEventRecord(stopEvent_inc,0);  //ending timing for inclusive
  	cudaEventSynchronize(stopEvent_inc);   
	cudaEventElapsedTime(&time, startEvent_inc, stopEvent_inc);   
 

	printf("%d\n%d\n", nV_1, nF_1);
	printf("center is %lf %lf %lf\n", center[0], center[1], center[2]);
	printf("normal of one direction is %lf %lf %lf\n", normal_1[0], normal_1[1], normal_1[2]);
	printf("%lf\n", vol);
	//free resources 
//	free(in); free(out); free(cuda_out);
	return 0;
}
