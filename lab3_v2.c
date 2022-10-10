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
	
    for (j=0; j<100; j++)  /*100 СЌРєСЃРїРµСЂРёРјРµРЅС‚РѕРІ */
    {
        int i = 0;
		unsigned int seed = j;
		
		double sum = 0;
		double del = 0;
		
		#pragma omp parallel default(none) private(i) shared(M1,M2,M3,N,seed,sum,del)
		{
       
			#pragma omp single
			for (i = 0; i<N; i++)
			{
				M1[i] = rand_r(&seed)%420;
			}
			
			#pragma omp single
			for (i = 0; i<N/2; i++)
			{
				double value = 420+rand_r(&seed)%3780;
				M3[i] = M2[i] = value;
			}
		
			#pragma omp for 
			for (i = 0; i<N; i++)
			{
				M1[i] = tanh(M1[i]);
			}
			
			#pragma omp barrier
			
			#pragma omp single
			M2[0] = log(fabs(tan(M2[0])));
			
			#pragma omp barrier
			
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
			
			#pragma omp barrier
			
			//Can't be parallelized cause it's an operation of sorting
			#pragma omp single
			heapSort(M2,N/2);
			
			#pragma omp barrier
			
			//Can't be parallelized cause it's an operation of searching FIRST not Null element
			#pragma omp single
			for (i = 0; i<N/2; i++) {
				if (M2[i]!=0) {
					del = M2[i];
					break;
				}
			}
			
			#pragma omp barrier
			
			#pragma omp for			
			for (i = 0; i<N/2; i++) {
				if ((int)(M2[i]/del)%2 == 0) {
					#pragma omp atomic
					sum += sin(M2[i]);
				}
			}
		}
		
       // printf("\nN=%d. %f", N, sum);
    }
    gettimeofday(&T2, NULL); 
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) /1000;
    printf("\nN=%d. Milliseconds passed: %ld", N, delta_ms); /* T2 - T1 */
    return 0;
}
