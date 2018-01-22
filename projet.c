#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

double* initializeGrid(int nbRow , int nbCol);
void drawGrid(double* grid, int nbRow , int nbCol);
double* computeGridSeq(double* grid, int nbRow, int nbCol, double dt, double dx);
double* computeGridPara(double* grid, int nbRow, int nbCol, double dt, double dx, int rang, int nbProcess, int nbRowPerProcess, int nbRowLastBlock);
int getNbRowPerProcess(int nbRow, int nbProcess, int rang);

/* Request 1D grid from 2D coordinates */
#define GRID(tab,i,j)  tab[j*nbCol+i]
// grid(tab, numéro de colonne , numéro de ligne)

/* Simplify the writing for cell's neighbours */
#define UP(tab,i,j) 	  ((j>0) ? GRID(tab,i,(j-1)) : 0)
#define DOWN(tab,i,j)     ((j<nbRow-1) ? GRID(tab,i,(j+1)) : 0)
#define RIGHT(tab,i,j)    ((i<nbCol-1) ? GRID(tab,(i+1),j) : 0)
#define LEFT(tab,i,j)     ((i>0) ? GRID(tab,(i-1),j) : 0)
// #define LEFTUP(tab,i,j)   ((i>0) && (j>0) ? GRID(tab,(i-1),(j-1)) : 0)
// #define RIGHTUP(tab,i,j)   ((i<nbCol-1) && (j>0) ? GRID(tab,(i+1),(j-1)) : 0)
// #define LEFTDOWN(tab,i,j)  ((i>0) && (j<nbRow-1) ? GRID(tab,(i-1),(j+1)) : 0)
// #define RIGHTDOWN(tab,i,j)  ((i<nbCol-1) && (j<nbRow-1) ? GRID(tab,(i+1),(j+1)) : 0)

static double WALLVALUE = -3.1415;

int main (int argc, char** argv){
  int rang, nbProcess;
  int nbRow = 30, nbCol = 20;
  double* globalGrid;
  double dt = 0.001, dx = 2*1/(double)nbCol;
  int speed = 100; // facteur par lequel diviser le temps de diffusion (ouai, si on augmente speed c'est plus lent, c'est logique)

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rang);
  MPI_Comm_size(MPI_COMM_WORLD, &nbProcess);

  if(rang == 0){
    printf("%d :: Initialize Grid\n",rang);

    globalGrid = initializeGrid(nbRow, nbCol);
    drawGrid(globalGrid, nbRow, nbCol);
  }

  int i;
  for(i = 0 ; i < 10000 ; i++){

    /*if(rang==0){
      globalGrid=computeGridSeq(globalGrid,nbRow,nbCol,dt,dx);
    }*/

    int nbRowPerProcess = getNbRowPerProcess(nbRow, nbProcess, rang);
    int nbRowLastBlock = getNbRowPerProcess(nbRow, nbProcess, nbProcess - 1);

    globalGrid = computeGridPara(globalGrid, nbRow, nbCol, dt, dx, rang, nbProcess, nbRowPerProcess, nbRowLastBlock);

    if(rang == 0){
      drawGrid(globalGrid, nbRow, nbCol);
      usleep(dt * pow(10,6) * speed);
    }
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}

double* initializeGrid(int nbRow , int nbCol){

  double heatPoint1[] = {8, 3, 8000};
  double heatPoint2[] = {18, 29, -100};

  double wallPoint1[]={5,5,WALLVALUE};
  double wallPoint2[]={6,5,WALLVALUE};
  double wallPoint3[]={7,5,WALLVALUE};
  double wallPoint5[]={8,5,WALLVALUE};
  double wallPoint6[]={9,5,WALLVALUE};
  double wallPoint7[]={10,5,WALLVALUE};
  double wallPoint8[]={11,5,WALLVALUE};
  double wallPoint4[]={12,5,WALLVALUE};



  double* grid = (double*) calloc(nbRow*nbCol, sizeof(double));

  GRID(grid, (int)heatPoint1[0], (int)heatPoint1[1]) = heatPoint1[2];
  GRID(grid, (int)heatPoint2[0], (int)heatPoint2[1]) = heatPoint2[2];
  GRID(grid, (int)wallPoint1[0], (int)wallPoint1[1]) = wallPoint1[2];
  GRID(grid, (int)wallPoint2[0], (int)wallPoint2[1]) = wallPoint2[2];
  GRID(grid, (int)wallPoint3[0], (int)wallPoint3[1]) = wallPoint3[2];
  GRID(grid, (int)wallPoint4[0], (int)wallPoint4[1]) = wallPoint4[2];
  GRID(grid, (int)wallPoint5[0], (int)wallPoint5[1]) = wallPoint5[2];
  GRID(grid, (int)wallPoint6[0], (int)wallPoint6[1]) = wallPoint6[2];
  GRID(grid, (int)wallPoint7[0], (int)wallPoint7[1]) = wallPoint7[2];
  GRID(grid, (int)wallPoint8[0], (int)wallPoint8[1]) = wallPoint8[2];

  return grid;
}

void drawGrid(double* grid, int nbRow , int nbCol){
  int col, row;

  system("clear");

  for(row = 0 ; row < nbRow ; row++){
    for(col = 0 ; col < nbCol ; col++){
      if((int)GRID(grid,col,row)==0){
        printf("      ");
      }else if(GRID(grid,col,row)==WALLVALUE){
        printf("++++++");
      }else{
        printf("%6.0f", GRID(grid, col, row));
      }
    }

    printf("\n");
  }
}

double* computeGridSeq(double* grid, int nbRow , int nbCol, double dt, double dx){
  double *newGrid = malloc(nbRow*nbCol*sizeof(double));
  int col, row;

  for(row = 0 ; row < nbRow ; row++){
    for(col = 0 ; col < nbCol ; col++){
      if(GRID(grid,col,row)!=WALLVALUE){
        double sumNeigbours = 0;
        int numberNeighbours=0;

        // Tests pour les cas aux limites de la grille
        if(row > 0) {
          if(UP(grid,col,row)!=WALLVALUE){
            sumNeigbours += UP(grid, col, row);
            numberNeighbours++;
          }
        }
        if(row < nbRow - 1) {
          if(DOWN(grid,col,row)!=WALLVALUE){
            sumNeigbours += DOWN(grid, col, row);
            numberNeighbours++;
          }
        }
        if(col > 0) {
          if(LEFT(grid,col,row)!=WALLVALUE){
            sumNeigbours += LEFT(grid, col, row);
            numberNeighbours++;
          }
        }
        if(col < nbCol - 1) {
          if(RIGHT(grid,col,row)!=WALLVALUE){
            sumNeigbours += RIGHT(grid, col, row);
            numberNeighbours++;
          }
        }

        //+ LEFTUP(grid,col,row) + RIGHTUP(grid,col,row) + LEFTDOWN(grid,col,row) + RIGHTDOWN(grid,col,row)
        GRID(newGrid, col, row) = GRID(grid, col, row) + (dt/dx/dx) * (sumNeigbours - numberNeighbours * GRID(grid, col,row));
      }else{
        GRID(newGrid, col, row) = GRID(grid, col, row);
      }
    }
  }

  return newGrid;
}

double* computeGridPara(double* grid, int nbRow , int nbCol, double dt, double dx, int rang, int nbProcess, int nbRowPerProcess, int nbRowLastBlock){

  int tag = 0;
  int process;
  double *localGrid = malloc((nbRowPerProcess+2) * nbCol*sizeof(double));
  double *newGrid = malloc(nbRow*nbCol*sizeof(double));

  // Scatter
  if(rang == 0){
    // Garde le premier bloc sur le process 0
    memcpy(localGrid, grid, (nbRowPerProcess + 1) * nbCol*sizeof(double)) ;

    // Envoie les blocs 1 à nbProcess-1
    for(process = 1 ; process <= nbProcess - 1 ; process++) {
      if(process != nbProcess - 1){
        MPI_Send(grid + (process * nbRowPerProcess * nbCol - nbCol), (nbRowPerProcess+2) * nbCol, MPI_DOUBLE, process, tag, MPI_COMM_WORLD);
      }
      else{
        MPI_Send(grid+(process*nbRowPerProcess*nbCol - nbCol), (nbRowLastBlock+1)*nbCol, MPI_DOUBLE, process, tag, MPI_COMM_WORLD);
      }
    }
  }
  else{
    if(rang != nbProcess - 1){
      MPI_Recv(localGrid, (nbRowPerProcess+2)*nbCol, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else{
      MPI_Recv(localGrid, (nbRowPerProcess+1)*nbCol, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  // Compute
  if(rang == 0 || rang == nbProcess - 1){
    localGrid = computeGridSeq(localGrid,nbRowPerProcess+1,nbCol,dt,dx);
  } else {
    localGrid = computeGridSeq(localGrid,nbRowPerProcess+2,nbCol,dt,dx);
  }

  // Gather
  double* localGridNoGhostZone = malloc(nbRowPerProcess * nbCol * sizeof(double));
  if(rang == 0){
    memcpy(newGrid, localGrid, nbRowPerProcess * nbCol * sizeof(double));

    for(process = 1 ; process <= nbProcess - 1 ; process++) {
      if(process != nbProcess - 1){
        MPI_Recv(newGrid + (process * nbRowPerProcess * nbCol), (nbRowPerProcess)*nbCol, MPI_DOUBLE, process, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      else {
        MPI_Recv(newGrid + (process * nbRowPerProcess * nbCol), (nbRowLastBlock)*nbCol, MPI_DOUBLE, process, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  }
  else{
    memcpy(localGridNoGhostZone, localGrid+nbCol, nbRowPerProcess*nbCol*sizeof(double));
    MPI_Send(localGridNoGhostZone, (nbRowPerProcess)*nbCol, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  }

  // MPI_Gather(localGridNoGhostZone, nbRowPerProcess*nbCol, MPI_DOUBLE, newGrid, nbRowPerProcess*nbCol, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return newGrid;
}

int getNbRowPerProcess(int nbRow, int nbProcess, int rang) {
  if(rang == nbProcess - 1) {
    return nbRow - (((int) nbRow / nbProcess) * (nbProcess - 1));
  }
  else {
    return (int) nbRow / nbProcess;
  }
}
