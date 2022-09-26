#include <stdio.h>
	#include <stdlib.h>
	#include <sys/time.h>
	#include <math.h>
	
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
	    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */
	    double *M1;  // указатель на массив
	    M1 = (double*)malloc(N * sizeof(double));
	    double *M2;  // указатель на массив
	    M2 = (double*)malloc(N/2 * sizeof(double));
	    double *M3;  // указатель на массив
	    M3 = (double*)malloc(N/2 * sizeof(double));
	    
	    int K;
	    K = atoi(argv[2]);
	    double *M1_temp;  // указатель на массив
	    M1_temp = (double*)malloc(N * sizeof(double));
	    fwSetNumThreads(K);

		
	    for (j=0; j<100; j++)  /*100 экспериментов */
	    {
	        srand(j);  /*инициализировать начальное значение ГСЧ */
			int i;
			unsigned int state;
	        for (i = 0; i<N; i++)
	        {
	            M1[i] = rand_r(&state)%420;
	        }
	        for (i = 0; i<N/2; i++)
	        {
	            double value = 420+rand_r(&state)%3780;
	            M2[i] = value;
		    if (i != 0) {
			M3[i] = M2[i-1];
			}
		    else {
		    	M3[i] = 0;
		    }
	        }
		
		fwsTanh_64f_A50(M1, M1_tmp, N);
	        //for (i = 0; i<N; i++)
	        //{
	        //    M1[i] = tanh(M1[i]);
	        //}
		
		fwsAdd_64f_I(M3, M2, N/2);
		fwsTan_64f_A50(M2, M2, N/2);
		fwsAbs_64f(M2, M2, N/2);
		fwsLn_64f_I(M2, N/2);
		fwsMul_64f_I(M1, M2, N/2);
	
	        //M2[0] = log(fabs(tan(M2[0])));
	        //for (i = 1; i<N/2; i++)
	        //{
	        //    M2[i] += M3[i-1];
	        //    M2[i] = log(fabs(tan(M2[i])));
	        //}
	        //for (i = 0; i<N/2; i++)
	        //{		
	        //    M2[i] = M2[i] * M1[i];
	        //}


	        heapSort(M2,N/2);
	        double del = 0;
	        for (i = 0; i<N/2; i++) {
	            if (M2[i]!=0) {
	                del = M2[i];
	                break;
	            }
	        }
	        double sum = 0;
	        for (i = 0; i<N/2; i++) {
	            if ((int)(M2[i]/del)%2 == 0) {
	                sum += sin(M2[i]);
	            }
	        }
	    }
	    gettimeofday(&T2, NULL); /* запомнить текущее время T2 */
	    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) /1000;
	    printf("\nN=%d. Milliseconds passed: %ld", N, delta_ms); /* T2 - T1 */
	    return 0;
	}

