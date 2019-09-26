#include <mpi.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include<math.h>


#define ROOT 0
#define BIG_NUMBER 10000

//Points struct:
struct
{
	int* xArr;
	int* yArr;
	int size;
}
typedef Points;


//-------------------Declaration of fuctions---------------------//
void initRootProcess(Points* allPoints, int size);
void fillWithNumbers(int* arrX, int* arrY, int size);
double handleMassive(Points* myPoints, int myPointsSize);
void checkAllocation(const void* p);
void handleLeftOvers(Points* allPoints, Points* myPoints, int numOfProcesses);
double massive(int x, int y);

//-----------------------------------------Main----------------------------------------//
int main(int argc, char *argv[])
{
	const int N = 32;
	int numOfProcs, myId, myPointsSize = 0, allPointsSize = 0;
	double t0 = 0, t1 = 0, mySum = 0, globalSum = 0;
	Points allPoints, myPoints;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &myId);

	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);

	if (myId == ROOT)
	{
		//initialize the array of points:
		initRootProcess(&allPoints, N);
		allPointsSize = allPoints.size;
		t0 = MPI_Wtime(); //start measuring time
	}

	//Every process gets the same number of points:
	myPointsSize = N * N / numOfProcs;

	//Root process will also be responsible for leftovers so it needs a bigger buffer:
	if (myId == ROOT)
		myPointsSize += N * N % numOfProcs;

	myPoints.xArr = (int*)calloc(myPointsSize, sizeof(int));
	checkAllocation(myPoints.xArr);

	myPoints.yArr = (int*)calloc(myPointsSize, sizeof(int));
	checkAllocation(myPoints.yArr);

	myPoints.size = myPointsSize;

	//Sends to each process X and Y values:
	MPI_Scatter(allPoints.xArr, N*N / numOfProcs, MPI_INT, myPoints.xArr, N*N / numOfProcs, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Scatter(allPoints.yArr, N*N / numOfProcs, MPI_INT, myPoints.yArr, N*N / numOfProcs, MPI_INT, ROOT, MPI_COMM_WORLD);

	//Add excess points to root:
	if (myId == ROOT)
		handleLeftOvers(&allPoints, &myPoints, numOfProcs);

	//every process calculates its points, and sums them up:
	mySum = handleMassive(&myPoints, myPointsSize);

	//sum all mySums to globalSum:
	MPI_Reduce(&mySum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);

	if (myId == ROOT)
	{
		t1 = MPI_Wtime(); //stop measuring time.

		printf("\nTotal time: %lf\n", t1 - t0);
		printf("Result is %e", globalSum);

		free(allPoints.xArr);
		free(allPoints.yArr);
	}

	free(myPoints.xArr);
	free(myPoints.yArr);

	MPI_Finalize();
}


//-------------------------------------Fuctions------------------------------------------//

double massive(int x, int y) {
	int i, loop = 1;
	double sum = 0;

	if (x < 5 || y < 5)
		loop = 10 + x + y;
	for (i = 0; i < loop* BIG_NUMBER + BIG_NUMBER % loop; i++)
		sum += cos(tan(sqrt(exp(sin((double)i / BIG_NUMBER)))));
	return sum;
}

void initRootProcess(Points* allPoints, int size)
{
	allPoints->xArr = (int*)calloc(size*size, sizeof(int));
	checkAllocation(allPoints->xArr);
	allPoints->yArr = (int*)calloc(size*size, sizeof(int));
	checkAllocation(allPoints->yArr);
	allPoints->size = size * size;
	fillWithNumbers(allPoints->xArr, allPoints->yArr, size);
}

void fillWithNumbers(int* arrX, int* arrY, int size)
{
	int i , j;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++) {
			arrX[j + i * size] = i;
			arrY[j + i * size] = j;
		}
	}
}

double handleMassive(Points* myPoints, int myPointsSize) {
	int i;
	double mySum = 0;
	for (i = 0; i < myPointsSize; i++)
	{
		mySum += massive(myPoints->xArr[i], myPoints->yArr[i]);
	}
	return mySum;
}

void handleLeftOvers(Points* allPoints, Points* myPoints, int numOfProcesses)
{
	int numOfLeftOvers = allPoints->size % numOfProcesses;
	int withoutLeftOvers = allPoints->size - numOfLeftOvers;
	memcpy(myPoints->xArr + allPoints->size / numOfProcesses,
		allPoints->xArr + withoutLeftOvers,
		numOfLeftOvers * sizeof(int));
	memcpy(myPoints->yArr + allPoints->size / numOfProcesses, allPoints->yArr + withoutLeftOvers,
		numOfLeftOvers * sizeof(int));
}

void checkAllocation(const void* p)
{
	if (!p){
		printf("Allocation failed.\n");
		exit(0);
	}
	return;
}