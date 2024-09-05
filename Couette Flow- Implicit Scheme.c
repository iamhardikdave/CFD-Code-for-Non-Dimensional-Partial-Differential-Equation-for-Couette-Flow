//Hardik Dave
// 234103324
//Impicit Scheme - BTCS TDMA

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h> // for RUN TIME Calculation
#include <unistd.h>
#define billion 1000000000.0
int main()
{
    struct timespec start, end;
    int n = 101;
    double dx = 1.0 / (n - 1);
    double Re = 100;
    double u[n], u_prev[n];
    double dt = 1e-2;
    double gamma = (dt / (Re * dx * dx));

    double err, a, b, c, d[n - 2], p[n - 2], q[n - 2];
    double x, k;
    FILE *fp1, *fp2, *fp3;
    fp1 = fopen("BTCS_TDMA time vs error.dat", "w");
    fp2 = fopen(" BTCS TDMA Log(error) vs Time.dat", "w");
    fp3 = fopen("All_data_file.txt", "w");

    int i;
    double t = 0.0;

    a = (1 + (2 * gamma));

    b = -gamma;

    c = -gamma;

    clock_gettime(CLOCK_REALTIME, &start);

    // Boundary conditions
    u[n - 1] = 1;
    u[0] = 0;
    for (i = 1; i < (n - 1); i++)
    {
        u[i] = 0; // Intial guess
    }

    // BTCS TDMA Formulation
    int iter = 0;
    char btcs_tdma[50];

    while (t < 50.005)
    {
        err = 0.0;
        if ((iter + 100) % 100 == 0)
        {
            sprintf(btcs_tdma, "vel-btcs_tdma_%.2lf.dat", t);
            FILE *fp4;
            fp4 = fopen(btcs_tdma, "w");
            for (i = 0; i < n; i++)
                fprintf(fp4, "%lf\t%lf\t\n", u[i], i * dx);
            fclose(fp4);
        }
        // Calculation of P & Q
        for (int i = 1; i < n - 1; i++)
        {
            if (i == 1)
            {
                d[i] = u[i];
                p[i] = -b / c;
                q[i] = d[i] / a;
            }
            else if (i < n - 2)
            {
                d[i] = u[i];
                p[i] = -b / (a + c * p[i - 1]);
                q[i] = (d[i] - c * q[i - 1]) / (a + c * p[i - 1]);
            }
            else if (i == n - 2)
            {
                d[i] = u[i] - (b * u[i + 1]);
                p[i] = 0;
                q[i] = (d[i] - c * q[i - 1]) / (a + c * p[i - 1]);
            }
        }
        // x sweep
        for (int i = n - 2; i > 0; i--)
        {
            k = u[i];
            if (i == n - 2)
            {
                u[i] = q[i];
            }
            else if (i < n - 2)
            {
                u[i] = p[i] * u[i + 1] + q[i];
            }
            err = err + pow((u[i] - k), 2);
        }

        err = sqrt((err / n));

		printf("iteration=%d\t", iter);
		printf("error=%lf\n", err);
		printf("time=%lf\t", t);
		printf("\n\n");

		fprintf(fp1, "%lf\t %lf\n", err, t);

               fprintf(fp2, "%lf \t %lf\n", log(err),t);

		fprintf(fp3, "iteration=%d \t error=%lf \n time=%lf\n\n", iter, err, t);
		iter++;
		t = t + dt;
	}
	fclose(fp1);
	fclose(fp2);


	clock_gettime(CLOCK_REALTIME, &end);
	// run_time = End - Start
	double run_time = (end.tv_sec - start.tv_sec) +
					  (end.tv_nsec - start.tv_nsec) / billion;

	printf("Run Time is %f seconds", run_time);
	fprintf(fp3, "Physical time taken is %f seconds", run_time);
	fclose(fp3);
}
