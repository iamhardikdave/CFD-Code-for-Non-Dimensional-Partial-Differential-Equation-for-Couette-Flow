//HARDIK DAVE
//234103324
// Couette flow problem using FTCS - Explicit Scheme

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>   // for RUN TIME Calculation
#include <unistd.h> 
#define billion 1000000000.0
int main()
{ 
  struct timespec start, end;
  
  int n = 101;   // No of grid points 
  int j ;
  double dy = 1.0 / (n - 1);
  double dt = 5e-3;  //time step
  double Re = 100;  //Reynolds no.
  double gamma = (dt / (Re * dy * dy));
  double u[n], u_prev[n];
  double t = 0.0;
  FILE *fp1, *fp2, *fp3, *fp4;
  fp1 = fopen("FTCS time vs error.dat", "w");
  fp2 = fopen("FTCS Log(error) vs Time.dat", "w");
  fp3 = fopen("All_data_file.txt", "w");
  clock_gettime(CLOCK_REALTIME, &start);
// Intialization and Boundary conditions
  for (j=0; j < n; j++)
  {
   if (j == 0)   //Bottom bc
   {
    u[j]= 0 ;
   }
   else if (j ==(n-1))   //Top bc
   {
    u[j]= 1 ;
   }
   else   //Interior points bc
   {
   u[j]= 0 ;
   }
  }

  float err;
  int iter = 0;
  char ftcs[50];

  while (t < 50.005)   //Assume steady state reached in 50 seconds
  {
    err = 0.0;

    if ((iter + 100) % 100 == 0)
    {
      sprintf(ftcs, "vel-ftcs_%.2f.dat", t);
      fp4 = fopen(ftcs, "w");
      for (j = 0; j < n; j++)
        fprintf(fp4, "%lf\t%lf\t\n", u[j], j * dy);
      fclose(fp4);
    }
    for (j = 0; j < n; j++)
    {
      u_prev[j] = u[j];
    }

    for (j = 1; j < (n - 1); j++)
    {
      u[j] = u_prev[j] + gamma * (u_prev[j + 1] - 2 * u_prev[j] + u_prev[j - 1]);
    }

    for (j = 0; j < n; j++)
      err = err + pow((u[j] - u_prev[j]), 2);

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
