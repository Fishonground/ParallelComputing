#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
//#include <omp.h>


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
            M2_v[i] = log(fabs(tan(M2_v[i] + M3_v[i - 1])));
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
    struct timeval T1, T2;
    long delta_ms;
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

    gettimeofday(&T1, nullptr); //запомнить текущее время T1

    //для сихронизации потоков
    cudaEvent_t syncEvent;
    cudaEventCreate(&syncEvent);    //Создаем event
    cudaEventRecord(syncEvent, nullptr);  //Записываем event

    //расчёт gridSize и blockSize для m1
    int gridSize, blockSize, minGridSize;
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, map_m1, 0, N);
    gridSize = (N + blockSize - 1) / blockSize;

    //расчёт gridSize и blockSize для m2
    int minGridSize2, blockSize2, gridSize2;
    cudaOccupancyMaxPotentialBlockSize(&minGridSize2, &blockSize2, map_m2, 0, N / 2);
    gridSize2 = (N / 2  + blockSize2 - 1) / blockSize2;
    
    
	
    for (j=0; j<100; j++)  /*100 экспериментов */
    {
        int i = 0;
		unsigned int seed = j;
		
		double sum = 0;
		double del = 0;
		
		//#pragma omp parallel default(none) private(i) shared(M1,M2,M3,N,seed,sum,del)
		//{
       
       
       		//GENERATE
       		
       		
			//#pragma omp single
			for (i = 0; i<N; i++)
			{
				M1[i] = rand_r(&seed)%420;
			}
			
			//#pragma omp single
			for (i = 0; i<N/2; i++)
			{
				M3[i] = M2[i] = 420+rand_r(&seed)%3780;
			}
			
			//копирование данных после инициализации
        		cudaMemcpy(M1_v, M1, sizeof(double) * N, cudaMemcpyHostToDevice);
        		cudaMemcpy(M2_v, M2, sizeof(double) * N / 2, cudaMemcpyHostToDevice);
        		cudaMemcpy(M3_v, M3, sizeof(double) * N / 2, cudaMemcpyHostToDevice);
			
		
			//MAP
			
			
			map_m1<<<gridSize, blockSize>>>(M1_v, N);
        		map_m2<<<gridSize2, blockSize2>>>(M2_v, M3_v, N / 2);
        		
        		cudaEventSynchronize(syncEvent);  //Синхронизируем event

        		cudaMemcpy(M1, M1_v, sizeof(double) * N, cudaMemcpyDeviceToHost);
        		cudaMemcpy(M2, M2_v, sizeof(double) * N / 2, cudaMemcpyDeviceToHost);
			
			//#pragma omp for 
			//for (i = 0; i<N; i++)
			//{
			//	M1[i] = tanh(M1[i]);
			//}
			
			//#pragma omp single
			M2[0] = log(fabs(tan(M2[0])));
			
			//#pragma omp for
			//for (i = 1; i<N/2; i++)
			//{
			//	M2[i] += M3[i-1];
			//	M2[i] = log(fabs(tan(M2[i])));
			//}
			
			//MERGE
			
			merge<<<gridSize2, blockSize2>>>(M1_v, M2_v, N / 2);

			cudaEventSynchronize(syncEvent);  //Синхронизируем event
			cudaMemcpy(M2, M2_v, sizeof(double) * N / 2, cudaMemcpyDeviceToHost);
			
			//#pragma omp for 
			//for (i = 0; i<N/2; i++)
			//{		
			//	M2[i] = M2[i] * M1[i];
			//}
			
			
			
			//SORT
			
			//Can't be parallelized cause it's an operation of sorting
			//#pragma omp single
			heapSort(M2,N/2);
		
			cudaEventSynchronize(syncEvent);  //Синхронизируем event
			
			//REDUCE
			
			
			//Can't be parallelized cause it's an operation of searching FIRST not Null element
			//#pragma omp single
			for (i = 0; i<N/2; i++) {
				if (M2[i]!=0) {
					del = M2[i];
					break;
				}
			}
			
			//#pragma omp for	reduction(+:sum)		
			for (i = 0; i<N/2; i++) {
				if ((int)(M2[i]/del)%2 == 0) {
					double x = sin(M2[i]);
					sum += x;
				}
			}
			
			cudaEventSynchronize(syncEvent);  //Синхронизируем event
		}
		
        //printf("\nN=%d. %f", N, sum);
    //}
    
    cudaEventDestroy(syncEvent);

    cudaFree(M1);
    cudaFree(M2);
    cudaFree(M3);

    free(M1);
    free(M2);
    free(M3);
    
    
    
    gettimeofday(&T2, NULL); /* запомнить текущее время T2 */
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) /1000;
    printf("\nN=%d. Milliseconds passed: %ld", N, delta_ms); /* T2 - T1 */
    return 0;
}
