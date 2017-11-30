// Reference Scan implementation - Author: Ananoymous student of ME759 Fall 2017
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include <cuda.h>
#include <math.h>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/scan.h>


void initializeArray(FILE* fp, float* ArrVer, float* ArrFace, int nVertices, int nFaces)
{
	char x[1024];
	long int tot_line = 0;
	while (fscanf(fp, "%1023s", x)) {
		tot_line++;
		if (tot_line > 20) break;
		puts(x);
        	if (x=="ply") continue;
		if ((x=="format")||(x=="comment")||(x=="property")) {fgets(x, 1023, fp);continue;}
								//	{fscanf(fp,"%[^\n]", x);continue;}
		if (x=="element") {
			fscanf(fp, "%1023s", x);
			if (x=="vertex") fscanf(fp, "%d", nVertices);
			else if (x=="face") fscanf(fp, "%d", nFaces);
		}
    	}
	
//	for( int i=0; i<nElements; i++){
//		int r=fscanf(fp,"%f",&arr[i]);
//		if(r == EOF){
//			rewind(fp);
//		}
//	}
}



int main(int argc, char* argv[]) {
	FILE *fp = fopen("bun_zipper.ply","r");
	//allocate resources
	int nV=0,nF=0;
	float* vertices, *faces;
	float time = 0.f;
	initializeArray(fp,vertices, faces,nV, nF);
	// Your code here
 
  
 // cudaMalloc((void**)&dout,size);
	//cudaMalloc((void**)&din,size);
 
  cudaEvent_t startEvent_inc, stopEvent_inc;
	cudaEventCreate(&startEvent_inc);
	cudaEventCreate(&stopEvent_inc);
  cudaEventRecord(startEvent_inc,0); // starting timing for inclusive  
  

	
  cudaEventRecord(stopEvent_inc,0);  //ending timing for inclusive
  cudaEventSynchronize(stopEvent_inc);   
	cudaEventElapsedTime(&time, startEvent_inc, stopEvent_inc);   
 

	printf("%d\n%d\n",nV,nF);


	//free resources 
//	free(in); free(out); free(cuda_out);
	return 0;
}
