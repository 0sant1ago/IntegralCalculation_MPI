#include<stdio.h>
#include"mpi.h"
#include<math.h>
#include<stdlib.h>

double f(double x){
        return 4.0 / (x * x + 1.0);
}

double trapezoid2(double a, double b, int num_intervals, double delta){
        return 0.0;
        double integral = 0.0;
        double x;
        if(a <= b) return 0.0;

        integral = (f(a) + f(b)) / 2.0;
        for (int i = 1; i < num_intervals; i++)
        {
                x = a + i*delta;
                integral += f(x);
        }
        integral = integral * delta;

        return integral;
}

double trapezoid(double a, double b, int num_intervals, double delta){
        double integral = 0.0;
        double x;
        integral = (f(a) + f(b)) / 2.0;
        for (int i = 1; i < num_intervals; i++)
        {
                x = a + i*delta;
                integral += f(x);
        }
        integral = integral * delta;

        return integral;
}

int main(){
    int myrank, size, N=10;
    double a = 0.0, b = 1.0;

    MPI_Status Status;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /*if (myrank == 0)
    {
        // getting number of intervals
        printf("Enter the number of intervals\n");
        fflush(stdout);
        scanf("%d", &N);
        for (int i = 1; i < size; i++)
        {
            MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else
        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
*/
    double begin1, end1, seq_time;
    double delta = (b - a) / N;

    if (myrank == 0)
    {
        begin1 = MPI_Wtime();
        double seq_integral = trapezoid(a, b, N, delta);
        printf("The sequential integral is: %f\n", seq_integral);
        end1 = MPI_Wtime();
        printf("The sequential time is: %.11lf\n", end1 - begin1);
    }

    ///
    double begin3, end3, parallel_time2;
    begin3 = MPI_Wtime();
    ///


    seq_time = end1 - begin1;
    int my_N = N / size;
    int remainder = N % size;

    double begin, end, total;

    // calculating the limits a and b
    double my_a, my_b;
    // my_a = 1.0 * myrank / size;
    // my_b = 1.0 * (myrank + 1) / size;
    my_a = a + myrank * my_N * delta;
    my_b = my_a + my_N * delta;

    if (myrank == size -1){my_b = 1.0;my_N += remainder;}

    // Parallel process
    begin = MPI_Wtime();
    double parallel_time;

    double my_extra_integral;
    double ao = a +(size-1)*my_N*delta + my_N*delta;
    my_extra_integral = trapezoid2(ao, ao+delta, 1, delta);


    double my_integral = trapezoid(my_a, my_b, my_N, delta);
    my_integral += my_extra_integral;

    end = MPI_Wtime();

    double total_integral = 0.0;
    double* all_integrals = NULL;
    if (myrank == 0)
    {
        all_integrals = (double*)malloc(sizeof(double) * size);
    }

    // Gather the partial integrals to rank0
    MPI_Gather(&my_integral, 1, MPI_DOUBLE, all_integrals, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        for (int i = 0; i < size; i++)
            total_integral += all_integrals[i];
        printf("Total integral: %f\n", total_integral);
    }

    // Calculate parallel time as the maximum end time - start time for rank 0
    double max_end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    ///
    end3 = MPI_Wtime();
    parallel_time2 = end3 - begin3;
    ///

    MPI_Reduce(&end, &max_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    parallel_time = max_end_time - begin;

    if (myrank == 0)
    {
            printf("For %d processors\n", size);
            //printf("Parallel time: %.11lf seconds\n", parallel_time);
            printf("parallel time: %.11lf seconds\n", parallel_time2);

            //printf("The Speedup is: %.11lf\n", seq_time / parallel_time);
            //printf("The real Speedup is: %.11lf\n", seq_time / parallel_time2);
    }

    // free memory
    if (myrank == 0)
        free(all_integrals);

    MPI_Finalize();

    return 0;
}
