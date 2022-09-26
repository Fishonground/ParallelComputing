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
    // Инициализируем наибольший элемент как корень
    int l = 2*i + 1; // левый = 2*i + 1
    int r = 2*i + 2; // правый = 2*i + 2
    // Если левый дочерний элемент больше корня
    if (l < n && arr[l] > arr[largest])
        largest = l;
    // Если правый дочерний элемент больше, чем самый большой элемент на данный момент
    if (r < n && arr[r] > arr[largest])
        largest = r;
    // Если самый большой элемент не корень
    if (largest != i)
    {
        swap(&arr[i], &arr[largest]);
    // Рекурсивно преобразуем в двоичную кучу затронутое поддерево
        heapify(arr, n, largest);
    }
}
void heapSort(double arr[], int n)
{
    // Построение кучи (перегруппируем массив)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);
    // Один за другим извлекаем элементы из кучи
    for (int i=n-1; i>=0; i--)
    {
        // Перемещаем текущий корень в конец
        swap(&arr[0], &arr[i]);
        // вызываем процедуру heapify на уменьшенной куче
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
	
    for (j=0; j<100; j++)  /*100 экспериментов */
    {
        srand(j);  /*инициализировать начальное значение ГСЧ */
		int i;
		unsigned int state;
        for (i = 0; i<N; i++)
        {
            M1[i] = rand_r(&state)%420;
			//printf("%d.%d, m1: %f\n", j, i, M1[i]);
        }
        for (i = 0; i<N/2; i++)
        {
            double value = 420+rand_r(&state)%3780;
            M3[i] = M2[i] = value;
			
			//printf("%d.%d, m2: %f\n", j, i, value);
        }
	
		
        for (i = 0; i<N; i++)
        {
            M1[i] = tanh(M1[i]);
            //M1[i] = cosh(M1[i])+1;
        }
		
        M2[0] = log(fabs(tan(M2[0])));
        for (i = 1; i<N/2; i++)
        {
			//printf("m2 before +m3: %f\n", M2[i]);
			//printf("m3 before +m3: %f\n", M3[i-1]);
            M2[i] += M3[i-1];
			
			//printf("m2 before log fabs tanh: %f\n", M2[i]);
			//printf("m2 tanh: %f\n", tan(M2[i]));
			//printf("m2 fabs tanh: %f\n", fabs(tan(M2[i])));
		//printf("m2 log fabs tanh: %f\n", log(fabs(tan(M2[i]))));
            M2[i] = log(fabs(tan(M2[i])));
			//printf("after +m3: %f\n", M2[i]);
        }
        for (i = 0; i<N/2; i++)
        {
			//printf("before mult: %f %f\n", M1[i], M2[i]);			
            M2[i] = M2[i] * M1[i];
			//printf("before sort: %f\n", M2[i]);
        }
        heapSort(M2,N/2);
        double del = 0;
        for (i = 0; i<N/2; i++) {
			
			//printf("after sort: %f\n", M2[i]);
			
            if (M2[i]!=0) {
                del = M2[i];
                break;
            }
        }
        double sum = 0;
        for (i = 0; i<N/2; i++) {
            if ((int)(M2[i]/del)%2 == 0) {
				//printf("before sin %f\n", M2[i]);
				
                sum += sin(M2[i]);
				//printf("sin %f\n", sin(M2[i]));
            }
        }
        //printf("\nN=%d. %f", N, sum);
    }
    gettimeofday(&T2, NULL); /* запомнить текущее время T2 */
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) /1000;
    printf("\nN=%d. Milliseconds passed: %ld", N, delta_ms); /* T2 - T1 */
    return 0;
}
