#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

void swap(double *x, double *y){
    double t = *x;
    *x = *y;
    *y = t;
}
void heapify(double arr[], int n, int i)
{
    int largest = i;
    int l = 2*i + 1; // Р»РµРІС‹Р№ = 2*i + 1
    int r = 2*i + 2; // РїСЂР°РІС‹Р№ = 2*i + 2
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
	int MAX_VALUE = 3780+420;
	
    int j, N, amount;
    struct timeval T1, T2;
    long delta_ms;
    N = atoi(argv[1]); /* N СЂР°РІРµРЅ РїРµСЂРІРѕРјСѓ РїР°СЂР°РјРµС‚СЂСѓ РєРѕРјР°РЅРґРЅРѕР№ СЃС‚СЂРѕРєРё */
    
    amount = atoi(argv[2]);
    omp_set_num_threads(amount);
    gettimeofday(&T1, NULL); /* Р·Р°РїРѕРјРЅРёС‚СЊ С‚РµРєСѓС‰РµРµ РІСЂРµРјСЏ T1 */
    double *M1;  // СѓРєР°Р·Р°С‚РµР»СЊ РЅР° РјР°СЃСЃРёРІ
    M1 = (double*)malloc(N * sizeof(double));
    double *M2;  // СѓРєР°Р·Р°С‚РµР»СЊ РЅР° РјР°СЃСЃРёРІ
    M2 = (double*)malloc(N/2 * sizeof(double));
    double *M3;  // СѓРєР°Р·Р°С‚РµР»СЊ РЅР° РјР°СЃСЃРёРІ
    M3 = (double*)malloc(N/2 * sizeof(double));
	
	int size_of_sub_array = (N/2)/amount;
	int size_of_last_array = size_of_sub_array+((N/2)-(size_of_sub_array*amount));
	
    for (j=0; j<100; j++)  /*100 СЌРєСЃРїРµСЂРёРјРµРЅС‚РѕРІ */
    {
        int i = 0;
		unsigned int seed = j;
		
		double sum = 0;
		double del = 0;
		
		#pragma omp parallel default(none) private(i) shared(M1,M2,M3,N,seed,sum,del,amount,size_of_sub_array,size_of_last_array,MAX_VALUE)
		{
			int current_thread_num = omp_get_thread_num();
       
			#pragma omp single
			for (i = 0; i<N; i++)
			{
				M1[i] = rand_r(&seed)%420;
			}
			
			#pragma omp single
			for (i = 0; i<N/2; i++)
			{
				M3[i] = M2[i] = 420+rand_r(&seed)%3780;
			}
		
			#pragma omp for 
			for (i = 0; i<N; i++)
			{
				M1[i] = tanh(M1[i]);
			}
			
			#pragma omp single
			M2[0] = log(fabs(tan(M2[0])));
			
			#pragma omp for
			for (i = 1; i<N/2; i++)
			{
				M2[i] += M3[i-1];
				M2[i] = log(fabs(tan(M2[i])));
			}
			
			#pragma omp for 
			for (i = 0; i<N/2; i++)
			{		
				M2[i] = M2[i] * M1[i];
			}
			
		
			int size; 
				
			if (current_thread_num<amount-1) {
				size = size_of_sub_array;
			} else {
				size = size_of_last_array;
			}
				
			heapSort((M2+current_thread_num*size_of_sub_array), size);

			
			#pragma omp barrier
			
			#pragma omp single
			{
				int *sub_array_current_pointers = (int*)malloc(amount * sizeof(int));
			
				for (i = 0; i < amount; i++) {
					sub_array_current_pointers[i] = 0;
				}
				
				double* sorted = (double*)malloc((N/2) * sizeof(double));
				
				for (i = 0; i<(N/2); i++) {
					int minimum = MAX_VALUE;
					int sub_array_with_minimum = 0;
					
					int j;
					
					for (j = 0; j < amount; j++) {
						int sub_array_offset = sub_array_current_pointers[j]; 
						
						int sub_array_size;
						
						if (j < amount-1) {
							sub_array_size = size_of_sub_array;
						} else {
							sub_array_size = size_of_last_array;	
						}						
						
						if (sub_array_offset < sub_array_size) {
							double challenger = *(M2+j*size_of_sub_array+sub_array_offset);
						
							if ( challenger < minimum ) {
								minimum = challenger;
								sub_array_with_minimum = j;
							}
						
						}
						
					}
					
					sorted[i] = minimum;
					sub_array_current_pointers[sub_array_with_minimum]++;
				}
				
				M2 = sorted;
				
			}

		
			
			//Can't be parallelized cause it's an operation of searching FIRST not Null element
			#pragma omp single
			for (i = 0; i<N/2; i++) {
				if (M2[i]!=0) {
					del = M2[i];
					break;
				}
			}
			
			#pragma omp for	reduction(+:sum)		
			for (i = 0; i<N/2; i++) {
				if ((int)(M2[i]/del)%2 == 0) {
					double x = sin(M2[i]);
					sum += x;
				}
			}
		}
		
        //printf("\nN=%d. %f", N, sum);
    }
    gettimeofday(&T2, NULL); /* Р·Р°РїРѕРјРЅРёС‚СЊ С‚РµРєСѓС‰РµРµ РІСЂРµРјСЏ T2 */
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) /1000;
    printf("\nN=%d. Milliseconds passed: %ld", N, delta_ms); /* T2 - T1 */
    return 0;
}
