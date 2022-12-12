#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>

__global__ void map_m1(double* M1_v, int size) {
    //линейный индекс потока
    unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;
    //сколько один поток выполняет
    unsigned int threadsNum = blockDim.x * gridDim.x;
    for (unsigned int i = id; i < size; i += threadsNum) {
        M1_v[i] = tanh(M1_v[i]);
    }
}

__global__ void map_m2(double* M2_v, double* M3_v, int size) {
    unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int threadsNum = blockDim.x * gridDim.x;
    
    for (unsigned int i = id; i < size; i += threadsNum) {
        if (i == 0) {
            M2_v[i] = log(fabs(tan(M2_v[i])));
        } else {
            M2_v[i] = log(fabs(tan(M2_v[i] + M3_v[i-1])));
        }
    }
}

__global__ void merge(const double* M1_v, double* M2_v, int size) {
    unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int threadsNum = blockDim.x * gridDim.x;
    for (unsigned int i = id; i < size; i += threadsNum) {
        M2_v[i] = (double) M1_v[i] * M2_v[i];
    }
}


void swap(double *x, double *y){
    double t = *x;
    *x = *y;
    *y = t;
}
void heapify(double arr[], int n, int i)
{
    int largest = i;
    int l = 2*i + 1; // левый = 2*i + 1
    int r = 2*i + 2; // правый = 2*i + 2
    if (l < n && arr[l] > arr[largest])
        largest = l;
    if (r < n && arr[r] > arr[largest])
        largest = r;
    if (largest != i)
    {
        swap(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

void heapSort(double arr[], int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);
        
    for (int i=n-1; i>=0; i--)
    {
        swap(&arr[0], &arr[i]);
        heapify(arr, i, 0);
    }
}

int main(int argc, char* argv[])
{
    int j, N;
    struct timeval T1, T2, T3, T4, T5;
    long delta_ms,delta_ms2,delta_ms3,delta_ms4;
    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    
    //omp_set_num_threads(amount);
    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */
    double *M1;  // указатель на массив
    M1 = (double*)malloc(N * sizeof(double));
    double *M2;  // указатель на массив
    M2 = (double*)malloc(N/2 * sizeof(double));
    double *M3;  // указатель на массив
    M3 = (double*)malloc(N/2 * sizeof(double));
    
    //lab 6
    double *M1_v, *M2_v, *M3_v;
    //выделение памяти на устройстве
    cudaMalloc(&M1_v, sizeof(double) * N);
    cudaMalloc(&M2_v, sizeof(double) * N / 2);
    cudaMalloc(&M3_v, sizeof(double) * N / 2);
    //gettimeofday(&T1, nullptr); //запомнить текущее время T1
    //printf("%ld, %ld\n", T1.tv_sec, T1.tv_usec);
    //printf("T1: %ld, %ld\n", T1.tv_sec, T1.tv_usec);
	
    //для сихронизации потоков
    cudaEvent_t syncEvent;
    cudaEventCreate(&syncEvent);    //Создаем event
    cudaEventRecord(syncEvent, NULL);  //Записываем event
    //расчёт gridSize и blockSize для m1
    int gridSize, blockSize, minGridSize;
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, map_m1, 0, N);
    gridSize = (N + blockSize - 1) / blockSize;
    //расчёт gridSize и blockSize для m2
    int minGridSize2, blockSize2, gridSize2;
    cudaOccupancyMaxPotentialBlockSize(&minGridSize2, &blockSize2, map_m2, 0, N / 2);
    gridSize2 = (N / 2  + blockSize2 - 1) / blockSize2;
    float ttime = 0;
	
    for (j=0; j<100; j++)  /*100 экспериментов */
    {
        int i = 0;
	unsigned int seed = j;
	double sum = 0;
	double del = 0;
       gettimeofday(&T1, nullptr);
       //GENERATE
	for (i = 0; i<N; i++)
	{
		M1[i] = rand_r(&seed)%420;
	}
	
	for (i = 0; i<N/2; i++)
	{
		M3[i] = M2[i] = 420+rand_r(&seed)%3780;
	}
	
	
	gettimeofday(&T2, nullptr);
	
	//копирование данных после инициализации
	cudaMemcpy(M1_v, M1, sizeof(double) * N, cudaMemcpyHostToDevice);
	cudaMemcpy(M2_v, M2, sizeof(double) * N / 2, cudaMemcpyHostToDevice);
	cudaMemcpy(M3_v, M3, sizeof(double) * N / 2, cudaMemcpyHostToDevice);
	
	
	cudaEvent_t start, stop;
	
	cudaEventCreate (&start);
	cudaEventCreate (&stop);
	cudaEventRecord (start, 0);
	//MAP
	map_m1<<<gridSize, blockSize>>>(M1_v, N);
	map_m2<<<gridSize2, blockSize2>>>(M2_v, M3_v, N / 2);
	cudaEventSynchronize(syncEvent);  //Синхронизируем event
	cudaMemcpy(M1, M1_v, sizeof(double) * N, cudaMemcpyDeviceToHost);
	cudaMemcpy(M2, M2_v, sizeof(double) * N / 2, cudaMemcpyDeviceToHost);
	
	
	
	//gettimeofday(&T3, NULL);
	
    	//printf("T4: %ld, %ld\n", T4.tv_sec, T4.tv_usec);
	//MERGE
	merge<<<gridSize2, blockSize2>>>(M1_v, M2_v, N / 2);
	cudaEventSynchronize(syncEvent);  //Синхронизируем event
	cudaMemcpy(M2, M2_v, sizeof(double) * N / 2, cudaMemcpyDeviceToHost);
	
	cudaEventRecord(stop,0);
	cudaEventSynchronize (stop);
	cudaEventElapsedTime (&ttime, start, stop);
	//printf("%.32f\n",ttime);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	
	
	gettimeofday(&T3, nullptr);
	//printf("T2: %ld, %ld\n", T2.tv_sec, T2.tv_usec);
	
	//SORT
	heapSort(M2,N/2);
	cudaEventSynchronize(syncEvent); 
	
	
	//REDUCE
	for (i = 0; i<N/2; i++) {
		if (M2[i]!=0) {
			del = M2[i];
			break;
		}
	}
	for (i = 0; i<N/2; i++) {
		if ((int)(M2[i]/del)%2 == 0) {
			double x = sin(M2[i]);
			sum += x;
		}
	}
	cudaEventSynchronize(syncEvent);  //Синхронизируем event
	gettimeofday(&T4, nullptr);
	
	//printf("N=%d. %f\n", N, sum);
	//printf("test: %d\n", T2.tv_sec-T1.tv_sec);
    	delta_ms2 = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) /1000;
    	//printf("test: %f\n", delta_ms2);
    	delta_ms3 = 1000*(T3.tv_sec - T2.tv_sec) + (T3.tv_usec - T2.tv_usec) /1000;
    	//delta_ms4 = 1000*(T4.tv_sec - T3.tv_sec) + (T4.tv_usec - T3.tv_usec) /1000;
    	//delta_ms3 = 0;
    	//delta_ms4 = 0.0;
    	//printf("%.16f\n", ttime);
    	printf("\nN=%d. Milliseconds passed: %ld, %ld \n", N, delta_ms2, delta_ms3);
    	//printf("\nN=%d. Milliseconds passed: %ld, %ld, %ld\n", N, delta_ms2, delta_ms3, delta_ms4); /* T2 - T1 */
    }
		
    
    cudaEventDestroy(syncEvent);
    cudaFree(M1);
    cudaFree(M2);
    cudaFree(M3);
    free(M1);
    free(M2);
    free(M3);
    
    gettimeofday(&T5, NULL); /* запомнить текущее время T2 */
    //delta_ms = 1000*(T5.tv_sec - T3.tv_sec) + (T5.tv_usec - T3.tv_usec) /1000;
    delta_ms = 1000*(T5.tv_sec - T1.tv_sec) + (T5.tv_usec - T1.tv_usec) /1000;
    //delta_ms2 = 1000*(T4.tv_sec - T3.tv_sec) + (T4.tv_usec - T3.tv_usec) /1000;
    //delta_ms3 = 1000*(T5.tv_sec - T4.tv_sec) + (T5.tv_usec - T4.tv_usec) /1000;
    printf("\nN=%d. Milliseconds passed: %ld \n", N, delta_ms); /* T2 - T1 */
    return 0;
}
