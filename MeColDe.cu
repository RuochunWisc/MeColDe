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

void initializeArr(FILE* fp, float* ArrVerX, float ArrVerY, float* ArrVerZ, float* ArrFace, int V, int F)
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
		fscanf(fp, "%f", &ArrFace[3*i]);
		fscanf(fp, "%f", &ArrFace[3*i+1]);
		fscanf(fp, "%f", &ArrFace[3*i+2]);
	}
}



int main(int argc, char* argv[]) {
	FILE *fp = fopen("bun_zipper.ply","r");
	//allocate resources
	int nV=0, nF=0;
	float time = 0.f;
	initializeSize(fp,&nV, &nF);
	float *vertices_x= (float *)malloc(sizeof(float)*nV);
	float *vertices_y= (float *)malloc(sizeof(float)*nV);
	float *vertices_z= (float *)malloc(sizeof(float)*nV);
        float *faces   = (float *)malloc(sizeof(float)*nF*3);
	initializeArr(fp,vertices_x,vertices_y,vertices_z,faces,nV,nF);
	// Your code here
 
  
 // cudaMalloc((void**)&dout,size);
	//cudaMalloc((void**)&din,size);
 
  cudaEvent_t startEvent_inc, stopEvent_inc;
	cudaEventCreate(&startEvent_inc);
	cudaEventCreate(&stopEvent_inc);
  cudaEventRecord(startEvent_inc,0); // starting timing for inclusive  
	//AABB construction
	//divergence theorem to locate the volume  

	
  cudaEventRecord(stopEvent_inc,0);  //ending timing for inclusive
  cudaEventSynchronize(stopEvent_inc);   
	cudaEventElapsedTime(&time, startEvent_inc, stopEvent_inc);   
 

	printf("%d\n%d\n", nV, nF);


	//free resources 
//	free(in); free(out); free(cuda_out);
	return 0;
}
