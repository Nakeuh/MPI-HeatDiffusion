#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

double* initializeGrid(int nbLine , int nbCol);
void drawGrid(double* grid, int nbLine , int nbCol);
double* computeGrid(double* grid, int nbLine , int nbCol, double dt, double dx);

/* Request 1D grid from 2D coordinates */
#define GRID(tab,i,j)  tab[j*nbCol+i]
// grid(tab, numéro de colonne , numéro de ligne)

/* Simplify the writing for cell's neighbours */
#define UP(tab,i,j) 	  ((j>0) ? GRID(tab,i,(j-1)) : 0)
#define DOWN(tab,i,j)     ((j<nbRow-1) ? GRID(tab,i,(j+1)) : 0)
#define RIGHT(tab,i,j)    ((i<nbCol-1) ? GRID(tab,(i+1),j) : 0)
#define LEFT(tab,i,j)     ((i>0) ? GRID(tab,(i-1),j) : 0)
#define LEFTUP(tab,i,j)   ((i>0) && (j>0) ? GRID(tab,(i-1),(j-1)) : 0)
#define RIGHTUP(tab,i,j)   ((i<nbCol-1) && (j>0) ? GRID(tab,(i+1),(j-1)) : 0)
#define LEFTDOWN(tab,i,j)  ((i>0) && (j<nbRow-1) ? GRID(tab,(i-1),(j+1)) : 0)
#define RIGHTDOWN(tab,i,j)  ((i<nbCol-1) && (j<nbRow-1) ? GRID(tab,(i+1),(j+1)) : 0)

int main (int argc, char** argv){
  int rang,nbProcess;
  int nbLine = 20, nbCol=20;
  double* globalGrid;
  double dt=0.001, dx=2*1/(double)nbCol;

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
    for(i=0;i<500;i++){
      printf("\n");
      printf("%d :: Compute Grid %d\n",rang,i);

      globalGrid=computeGrid(globalGrid,nbLine,nbCol,dt,dx);
      system("clear");
      drawGrid(globalGrid,nbLine,nbCol);
      usleep(10000*10);
    }
  }
  MPI_Finalize();

  return EXIT_SUCCESS;
}

double* initializeGrid(int nbRow , int nbCol){
  // X, Y , Temperature
  //int heatPoint1[] = {3,0,50};

  // Interesting with a 3x3 grid
/*  int heatPoint1[] = {0,0,100};
  int heatPoint2[] = {1,0,200};
  int heatPoint3[] = {2,0,300};
  int heatPoint4[] = {0,1,400};
  int heatPoint5[] = {1,1,500};
  int heatPoint6[] = {2,1,600};
  int heatPoint7[] = {0,2,700};
  int heatPoint8[] = {1,2,800};
  int heatPoint9[] = {2,2,900};*/

  int heatPoint1[] = {1,2,8000};
  int heatPoint2[] = {18,19,15000};

  double* grid = (double*)calloc(nbRow*nbCol,sizeof(double));

  GRID(grid,heatPoint1[0],heatPoint1[1])= heatPoint1[2];
  GRID(grid,heatPoint2[0],heatPoint2[1])= heatPoint2[2];
  /*GRID(grid,heatPoint3[0],heatPoint3[1])= heatPoint3[2];
  GRID(grid,heatPoint4[0],heatPoint4[1])= heatPoint4[2];
  GRID(grid,heatPoint5[0],heatPoint5[1])= heatPoint5[2];
  GRID(grid,heatPoint6[0],heatPoint6[1])= heatPoint6[2];
  GRID(grid,heatPoint7[0],heatPoint7[1])= heatPoint7[2];
  GRID(grid,heatPoint8[0],heatPoint8[1])= heatPoint8[2];
  GRID(grid,heatPoint9[0],heatPoint9[1])= heatPoint9[2];
*/
  return grid;
}

void drawGrid(double* grid, int nbRow , int nbCol){
  int col,row;
  for(row=0; row<nbRow;row++){
    for(col=0; col<nbCol; col++){
      printf("%6.0f", GRID(grid,col,row));
    }
    printf("\n");
  }
}

double* computeGrid(double* grid, int nbRow , int nbCol, double dt, double dx){
  double *tmp = malloc(nbRow*nbCol*sizeof(double));
  int col,row;
  for(row=0; row<nbRow;row++){
    printf("Line %d\n",row);

    for(col=0; col<nbCol; col++){
      double sumNeigbours = 0;
      int numberNeighbours=0;

      // Tests pour les cas aux limites de la grille
      if(row>0){
        sumNeigbours+=UP(grid,col,row);
        numberNeighbours++;
      }
      if(row<nbRow-1){
        sumNeigbours+=DOWN(grid,col,row);
        numberNeighbours++;
      }
      if(col>0){
        sumNeigbours+=LEFT(grid,col,row);
        numberNeighbours++;
      }
      if(col<nbCol-1){
        sumNeigbours+=RIGHT(grid,col,row);
        numberNeighbours++;
      }

      //+ LEFTUP(grid,col,row) + RIGHTUP(grid,col,row) + LEFTDOWN(grid,col,row) + RIGHTDOWN(grid,col,row)
      GRID(tmp,col,row) = GRID(grid,col,row) + (dt/dx/dx)* (sumNeigbours-numberNeighbours*GRID(grid,col,row));
    }
  }

  return tmp;
}
