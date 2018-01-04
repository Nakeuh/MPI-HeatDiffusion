#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

double* initializeGrid(int nbLine , int nbCol);
void drawGrid(double* grid, int nbLine , int nbCol);
double* computeGrid(double* grid, int nbLine , int nbCol, double dt, double dx);

/* Request 1D grid from 2D coordinates */
#define GRID(tab,i,j)  tab[j*nbCol+i]


/* Simplify the writing for cell's neighbours */
#define UP(tab,i,j) 	  ((j>0) ? GRID(tab,i,j-1) : 0)
#define DOWN(tab,i,j)     ((j<nbLine-1) ? GRID(tab,i,j+1) : 0)
#define RIGHT(tab,i,j)    ((i<nbCol-1) ? GRID(tab,i+1,j) : 0)
#define LEFT(tab,i,j)     ((i>0) ? GRID(tab,i-1,j) : 0)
#define LEFTUP(tab,i,j)   ((i>0) && (j>0) ? GRID(tab,i-1,j-1) : 0)
#define RIGHTUP(tab,i,j)   ((i<nbCol-1) && (j>0) ? GRID(tab,i+1,j-1) : 0)
#define LEFTDOWN(tab,i,j)  ((i>0) && (j<nbLine-1) ? GRID(tab,i-1,j+1) : 0)
#define RIGHTDOWN(tab,i,j)  ((i<nbCol-1) && (j<nbLine-1) ? GRID(tab,i+1,j+1) : 0)

int main (int argc, char** argv){
  int rang,nbProcess;
  int nbLine = 20, nbCol=20;
  double* globalGrid;
  double dt=1, dx=1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rang);
  MPI_Comm_size( MPI_COMM_WORLD, &nbProcess);

  //int nbLinePerProcess = nbLine/nbProcess;

  if(rang==0){
    printf("%d :: Initialize Grid\n",rang);
    //printf("%d :: Nb line per process %d\n",rang, nbLinePerProcess);

    globalGrid=initializeGrid(nbLine,nbCol);
    drawGrid(globalGrid,nbLine,nbCol);

    int i;
    for(i=0;i<2;i++){
      printf("\n");
      printf("%d :: Compute Grid %d\n",rang,i);
      globalGrid=computeGrid(globalGrid,nbLine,nbCol,dt,dx);

      drawGrid(globalGrid,nbLine,nbCol);
    }
  }



  MPI_Finalize();

  return EXIT_SUCCESS;
}

double* initializeGrid(int nbLine , int nbCol){
  // X, Y , Temperature
  int heatPoint1[] = {0,0,50};
  int heatPoint2[] = {10,3,90};

  double* grid = calloc(nbLine*nbCol,sizeof(double));

  GRID(grid,heatPoint1[0],heatPoint1[1])= heatPoint1[2];
  GRID(grid,heatPoint2[0],heatPoint2[1])= heatPoint2[2];

  return grid;
}

void drawGrid(double* grid, int nbLine , int nbCol){
  int i,j;
  for(i=0; i<nbLine;i++){
    for(j=0; j<nbCol; j++){
      printf("%6.0f", GRID(grid,j,i));
    }
    printf("\n");
  }
}

double* computeGrid(double* grid, int nbLine , int nbCol, double dt, double dx){
  double *tmp = malloc(nbLine*nbCol*sizeof(double));
  int i,j;
  for(i=0; i<nbLine;i++){
    for(j=0; j<nbCol; j++){
      double sumNeigbours = UP(grid,j,i) + DOWN(grid,j,i) + LEFT(grid,j,i) + RIGHT(grid,j,i) + LEFTUP(grid,j,i) + RIGHTUP(grid,j,i) + LEFTDOWN(grid,j,i) + RIGHTDOWN(grid,j,i);
      GRID(tmp,j,i) = GRID(grid,j,i) + (dt/dx/dx)* (sumNeigbours-8*GRID(grid,j,i));
    }
  }

  return tmp;
}
