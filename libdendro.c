/**
 * The following procedures are C implementations of some loop-intensive dendrogram code. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "idl_export.h"
#define TESTING 0
#define MAX(a, b) (( (a) > (b)) ? (a) : (b))
#define MIN(a, b) (( (a) < (b)) ? (a) : (b))
#define OFF2(x, y, xsize) ((x) + ((y) * (xsize)))
#define OFF3(x, y, z, xsize, ysize) ((x) + ((y) * (xsize)) + ((z) * (xsize) * (ysize)))
#define STACK2(di,dj) off=OFF2(i+(di),j+(dj),xsize);			\
  if( (i+(di) >= 1) && (i+(di) < xsize - 1) &&				\
      (j+(dj) >= 1) && (j+(dj) < ysize - 1) &&				\
      *(result + off) == 0 && *(data + off) > thresh) {			\
    *(stack + ++top) = off;						\
    *(result + off) = 1;						\
  } else if(*(result + off) == 0) *(result + off) = 2

#define STACK3(di,dj,dk) off=OFF3(i+(di),j+(dj),k+(dk),xsize,ysize);	\
  if( (i+(di) >= 1) && (i+(di) < xsize - 1) &&				\
      (j+(dj) >= 1) && (j+(dj) < ysize - 1) &&				\
      (k+(dk) >= 1) && (k+(dk) < zsize - 1) &&				\
      *(result + off) == 0 && *(data + off) > thresh) {			\
    *(stack + ++top) = off;						\
    *(result + off) = 1;						\
  } else if(*(result + off) == 0) *(result + off) = 2

/**
 * Creates a local maxima mask in a 2d image.
 * Local maxima are greater than all of their neighbors in 
 * a square of size 2 * nfriend + 1
 *
 * Argv is a 5 element array:
 *   argv[0]: the data array, double precision
 *   argv[1]: nfriend, an IDL_LONG
 *   argv[2]: x size of the data array, an IDL_LONG
 *   argv[3]: y size of the data array, an IDL_LONG
 *   argv[4]: The output mask array. Type UCHAR (IDL byte). Assumed to be populated with 1s
 *
 * 
 */
void alllocmax_2d(IDL_INT argc, void *argv[]) {
  
  double *cube = argv[0];
  IDL_LONG friends = *(IDL_LONG *) argv[1];
  IDL_LONG xsize = *(IDL_LONG *) argv[2];
  IDL_LONG ysize = *(IDL_LONG *) argv[3];
  UCHAR *result = argv[4];

  IDL_LONG i,j,fail;
  IDL_LONG di, dj;
  IDL_LONG l, r, b, t;
  // loop over pixels
  for(j = 0; j < ysize; j++) {
    for(i = 0; i < xsize; i++) {
      int off = OFF2(i,j,xsize);
      l = MAX(0, i - friends);
      r = MIN(xsize - 1, i + friends);
      b = MAX(0, j - friends);
      t = MIN(ysize - 1, j + friends);
      fail = 0;
      // loop over neighbors
      for(di = l; di <= r; di++) {
	for(dj = b; dj <= t; dj++) {
	  if (di == i && dj == j) continue;
	  fail = (*(cube + off) <= *(cube + OFF2(di, dj, xsize)));
	  if (fail == 1) break;
 	}  // neighbor y loop
	if (fail == 1) break;
      } // neighbor x loop
      *(result + off) = fail ? 0 : 1;
    } // x loop
  } // y loop
}


/** 
 * Finds local maxima in a 3D array. l
 * Local maxima are pixels which are greater than all neighbors 
 * in a rectangle of radius (friends, friends, specfriends).
 *
 * argv is a 7 element array:
 *  argv[0]: Data array, double
 *  argv[1]: friends, IDL_LONG
 *  argv[2]: specfriends, IDL_LONG
 *  argv[3]: xsize, IDL_LONG
 *  argv[4]: ysize, IDL_LONG
 *  argv[5]: zsize, IDL_LONG
 *  argv[6]: result array, UCHAR (IDL Byte)
 */
void alllocmax_3d(IDL_INT argc, void *argv[]) {
  
  double *cube = argv[0];
  IDL_LONG friends = *(IDL_LONG *) argv[1];
  IDL_LONG specfriends = *(IDL_LONG *) argv[2];
  IDL_LONG xsize = *(IDL_LONG *) argv[3];
  IDL_LONG ysize = *(IDL_LONG *) argv[4];
  IDL_LONG zsize = *(IDL_LONG *) argv[5];
  UCHAR *result = argv[6];

  IDL_LONG i,j,k,fail;
  IDL_LONG di, dj, dk;
  IDL_LONG l, r, b, t, fr, ba;
  // loop over pixels
  for(k = 0; k < zsize; k++) {
    for(j = 0; j < ysize; j++) {
      for(i = 0; i < xsize; i++) {
	int off = OFF3(i,j,k,xsize,ysize);
	l = MAX(0, i - friends);
	r = MIN(xsize - 1, i + friends);
	b = MAX(0, j - friends);
	t = MIN(ysize - 1, j + friends);
	fr = MAX(0, k - specfriends);
	ba = MIN(zsize - 1, k + specfriends);
      fail = 0;
      // loop over neighbors
      for(di = l; di <= r; di++) {
	for(dj = b; dj <= t; dj++) {
	  for(dk = fr; dk <= ba; dk++) {
	    if (di == i && dj == j && dk == k) continue;
	    fail = (*(cube + off) <= *(cube + OFF3(di,dj,dk,xsize,ysize)));
	    if (fail == 1) break;
	  } // z neighbor
	  if (fail == 1) break;
	} // y neighbor
	if (fail == 1) break;
      } // x neighbor
      *(result + off) = fail ? 0 : 1;
      } // x pixel
    } // y pixel
  } // z pixel
} // alllocmax_3d

 
void fill_2d(double *data, IDL_LONG xsize, IDL_LONG ysize,
	     int i, int j, int all,
	     double thresh, UCHAR *result) {
  int off;        
  int top, capacity = xsize;
  off = OFF2(i,j,xsize);
  
  int *stack;

  // easy case: seed position is zero. quit
  if (*(data + off) <= thresh ||
      (i < 1) || (i >= xsize - 1) || 
      (j < 1) || (j >= ysize - 1)) {
    printf("XXX Quitting\n");
    return;
  }

  stack = malloc(capacity * sizeof(int));
  top = -1;

  *(stack + ++top) = OFF2(i,j, xsize);
  while(top >= 0) {
    // grow stack, if needed
    if(top + 10 > capacity) {
      int newcap = MIN(capacity * 2, xsize * ysize);
      //printf("growing %i\n", capacity);
      int *newstack;
      if(!(newstack = malloc(newcap * sizeof(int)))) {
	printf("Memory allocation failure in fill_2d. Aborting");
	return;
      }
      memcpy(newstack, stack, capacity * sizeof(int));
      free(stack);
      stack = newstack;
      capacity = newcap;
    }

    off = *(stack + top--);
    i = off % xsize; j = off / xsize;  
    
    *(result + off) = 1;
    
    // add neighbors to the stack
    STACK2(1,0); STACK2(-1,0); 
    STACK2(0,1); STACK2(0,-1);
    if (all != 0) {
      STACK2(1,1); STACK2(1,-1); 
      STACK2(-1,1); STACK2(-1,-1);
    }
  } // top >= 0
free(stack);
}

void fill_3d(double *data, 
	     IDL_LONG xsize, IDL_LONG ysize, IDL_LONG zsize, 
	     int i, int j, int k, int all, 
	     double thresh, UCHAR *result) {
  int off, top;
  off = OFF3(i,j,k,xsize,ysize);
  int *stack, capacity = xsize;
  // easy case: seed position is zero. quit.
  if ((i < 1) || (i >= xsize-1) || 
      (j < 1) || (j >= ysize-1) || 
      (k < 1) || (k >= zsize-1) ||
      *(data + off) <= thresh) return;
 
  stack = malloc(capacity * sizeof(int));
  if (stack == NULL) {
    printf("allocation failure\n");
    return;
  }
  top = -1;

  *(stack + ++top) = OFF3(i,j,k,xsize,ysize);
  while(top >= 0) {
    // grow the stack, if needed
    if(top + 30 > capacity) {
      int newcap = capacity * 2;
      int *newstack = malloc(newcap * sizeof(int));
      memcpy(newstack, stack, capacity * sizeof(int));
      free(stack);
      stack = newstack;
      capacity = newcap;
    }

    off = *(stack + top--);
    i = off % xsize; j = (off / xsize) % ysize; k = off / (xsize * ysize);
    *(result + off) = 1;
    
    STACK3(0,0,1); STACK3(0,0,-1);
    STACK3(0,1,0); STACK3(0,-1,0);
    STACK3(1,0,0); STACK3(-1,0,0);
    if (all != 0) {
      STACK3(0,1,1); STACK3(0,1,-1); STACK3(0,-1,1); STACK3(0,-1,-1);
      STACK3(1,0,1); STACK3(1,0,-1); STACK3(-1,0,1); STACK3(-1,0,-1);
      STACK3(1,1,0); STACK3(1,-1,0); STACK3(-1,1,0); STACK3(-1,-1,0);
      STACK3(1,1,1); STACK3(1,1,-1); STACK3(1,-1,1); STACK3(1,-1,-1);
      STACK3(-1,1,1); STACK3(-1,1,-1); STACK3(-1,-1,1); STACK3(-1,-1,-1);
    }
  }
  free(stack);
  return;
}
 
void fill(int argc, void* argv[]) {
  double *cube   = argv[0];
  UCHAR *result  = argv[1];
  IDL_LONG ndim  = *(IDL_LONG *) argv[2];
  IDL_LONG xsize = *(IDL_LONG *) argv[3];
  IDL_LONG ysize = *(IDL_LONG *) argv[4];
  IDL_LONG zsize = *(IDL_LONG *) argv[5];
  IDL_LONG xseed = *(IDL_LONG *) argv[6];
  IDL_LONG yseed = *(IDL_LONG *) argv[7];
  IDL_LONG zseed = *(IDL_LONG *) argv[8];
  double thresh  = *(double *) argv[9];
  IDL_LONG all   = *(IDL_LONG *) argv[10];
  int i;
  
  if (ndim == 2) {
    fill_2d(cube, xsize, ysize,(int)xseed, (int)yseed, (all == 0) ? 0 : 1, 
	    thresh, result);
    for(i = 0; i < xsize * ysize; i++) if (result[i] == 2) result[i] = 0;
  } else {
    fill_3d(cube, xsize, ysize, zsize, (int)xseed, (int)yseed, (int)zseed, 
	    (all == 0) ? 0 : 1, thresh, result);
    for(i = 0; i < xsize * ysize * zsize; i++) if (result[i] == 2) result[i] = 0;
  }
}
  
int main(int argc, char** argv) {
  int i, j;
  double data[25], thresh;
  UCHAR result[25];
  thresh = 1;

  for(i = 0; i < 25; i++) *(data + i) = 0;
  

  for (i = 1; i < 3; i++) {
    for (j = 1; j < 3; j++) {
      int off = OFF2(i,j, 5);
      *(data + off) = 3;
    }
  }

  for(i = 0; i < 25; i++) *(result + i) = 0;  

  fill_2d(data, 5, 5, 2,2,1,thresh,result);

  for(i = 0; i < 25; i++) printf("%f ", *(data + i));
  printf("\n");
  for(i = 0; i < 25; i++) printf("%i ", *(result + i));
  printf("\n");
  return 0;
}

