#include <stdio.h>
#include <sys/types.h>
#include <time.h>

extern "C"
{
  double tmrc(void)
  {
    int rc;
    //struct timeval tv;
    struct timespec tv;

    /*rc = gettimeofday (&tv, NULL);
    if (rc == -1) {
      fprintf (stderr,"error in tmrc for gettimeofday\n");
      return 0.;
    }*/

    clock_gettime(CLOCK_MONOTONIC, &tv);
    //return ((double) tv.tv_sec) + 1.e-6 * ((double) tv.tv_usec);

    return ((double) tv.tv_sec) + 1.e-9 * ((double) tv.tv_nsec);

    /* return time(NULL); */
    /* return clock()/CLOCKS_PER_SEC; */
  }
}
