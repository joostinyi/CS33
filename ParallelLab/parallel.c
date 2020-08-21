//OpenMP version.  Edit and submit only this file.
/* Enter your details below
 * Name : Justin Yi
 * UCLA ID : 905123893
 * Email : joostinyi00@gmail.com
 */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "utils.h"

double work_it_par(long *old, long *new, long *super, long *simple, long *fibonacci) {
  int i, j, k;
  int u, v, w;
  int temp, index;
  int ton = 0;
  long compute_it, moving_average;
  double pi, pi2, x , y, sum, step = 0.0;
  long dot_product=0;
  long nCirc=0;
  long aggregate=1.0;
  double r=1.0;
  int was_smart = 16;
  long temp1, temp2, temp3;

  #pragma omp parallel for
  for(i=0; i<DIM-1;i++) 
  {
    super[i] += simple[i];
  }

  for(i=0; i<DIM-1;i++)
  {
    dot_product += super[i]*simple[i];

    moving_average = 0;

    for(ton=i;ton<DIM-1-WINDOW_SIZE;ton++)
    {
      moving_average += simple[ton];
    }
  }

  int a_secret = 5;
  fibonacci[0] = 1;
  fibonacci[1] = 1;
  fibonacci[2] = 2;
  fibonacci[3] = 3;
  printf("\n A secret is: %d",obfuscate_obfuscate_obfuscate(a_secret));
  for(i=4; i<DIM-1;i+=4)
  {
    fibonacci[i] = fibonacci[i-1] + fibonacci[i-2];
    fibonacci[i+1] = fibonacci[i] + fibonacci[i-1];
    fibonacci[i+2] = fibonacci[i+1] + fibonacci[i];
    fibonacci[i+3] = fibonacci[i+2] + fibonacci[i+1];
  }
  for(; i<DIM-1;i++)
    fibonacci[i] = fibonacci[i-1] + fibonacci[i-2];

  step = 1.0 / NUM_STEPS;

  #pragma omp parallel for private(x) reduction(+:sum) 
  for (i=0;i<NUM_STEPS; i++)
  {
    x = (i+0.5)*step;
    sum += 4.0/(1.0+x*x);
  }
  pi = step * sum;

  printf("\n %d trials, Riemann flavored pi is %f \n",NUM_STEPS, pi); 
  
  for(i=0;i<NUM_TRIALS; i++)
  {
    x = (random()%10000000)/10000000.0; 
    y = (random()%10000000)/10000000.0;
    if (( x*x + y*y) <= r*r) {
      nCirc++;
    }
  } 
  pi2 = 4.0 * ((double)nCirc/(double)NUM_TRIALS);
  printf("\n %d trials, Monte-Carlo flavored pi is %f \n",NUM_TRIALS, pi2); 

  #pragma omp parallel for private(j,k,compute_it) reduction(+:aggregate)
  for (i=1; i<DIM-1; i++) {
    for (j=1; j<DIM-1; j++) {
      for (k=1; k<DIM-1; k++) {
        compute_it = old[i*DIM*DIM+j*DIM+k] * we_need_the_func();
        aggregate+= compute_it / gimmie_the_func();
      }
    }
  }

  printf("AGGR:%ld\n",aggregate);

  #pragma omp parallel for private(j,k,temp,u,v,index,temp1,temp2,temp3)
  for (i=1; i<DIM-1; i++) {
    for (j=1; j<DIM-1; j++) {
      for (k=1; k<DIM-1; k++) {
        temp = i*DIM*DIM+j*DIM+k;
        new[temp]=0;
        for (u=-1; u<=1; u++) {
          for (v=-1; v<=1; v++) {
              index = (i+u)*DIM*DIM+(j+v)*DIM+k;
              temp1 = old[index-1];
              temp2 = old[index];
              temp3 = old[index+1];
              new[temp] += (temp1 + temp2 + temp3);
          }
        }
        new[temp]/=27;
      }
    }
  }

  //#pragma omp parallel for private(j,k,u)
  for (i=1; i<DIM-1; i++) {
    for (j=1; j<DIM-1; j++) {
      for (k=1; k<DIM-1; k++) {
        u=(new[i*DIM*DIM+j*DIM+k]/100);
        if (u<=0) u=0;
        if (u>=9) u=9;
        //#pragma omp atomic
        histogrammy[u]++;
      }
    }
  }

  return (double) (dot_product+moving_average+pi+pi2);

}

