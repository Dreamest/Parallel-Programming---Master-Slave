#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//defines:
#define ROOT 0
#define DATA_TAG 0
#define TERMINATION_TAG 1
#define BIG_NUMBER 10000
#define TASK_SIZE 2

//Points struct:
struct
{
	int* xArr;
	int* yArr;
	int size;
}
typedef Points;

//-------------------Declaration of fuctions---------------------//
void masterProcess(int numOfProcesses);
void slaveProcess();
void initRootProcess(Points* allPoints, int size);
void fillWithNumbers(int* arrX, int* arrY, int size);
void checkAllocation(const void* p);
double massive(int x, int y);
int getNextTask(int* nextTask, Points* allPoints, int offset);

//------------------------Main--------------------------//
int main(int argc, char *argv[])
{
	int numProcs, myId;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	if (myId == ROOT)
	{
		masterProcess(numProcs);
	}
	else
	{
		slaveProcess();
	}
	MPI_Finalize();
}

//------------------------Master--------------------------//
void masterProcess(int numOfProcesses)
{
	double t0, t1;
	const int N = 32;
	MPI_Status status;
	int id, offset = 0;
	int jobsSent = 0, jobsToBeDone = N * N / (TASK_SIZE/2);
	double globalSum = 0;
	Points allPoints;
	int task[TASK_SIZE]; //task will be filled with an X value followed by a Y value

	//loads X and Y values to allPoints struct.
	initRootProcess(&allPoints, N);

	t0 = MPI_Wtime(); //start measuring time

	//Start by sending one task for each slave
	for (id = 1; id < numOfProcesses; id++)
	{
		offset = getNextTask(task, &allPoints, offset); //loads into task X and Y values based on offset, and update offset.
		MPI_Send(task, TASK_SIZE, MPI_INT, id, DATA_TAG, MPI_COMM_WORLD);
	}

	//In a loop, recieve result then send to the process that returned it a new task until all tasks are done:
	for (jobsSent = numOfProcesses - 1; jobsSent < jobsToBeDone; jobsSent++)
	{
		int source;
		double localSum;
		int tag = jobsToBeDone - jobsSent > numOfProcesses - 1 ? DATA_TAG : TERMINATION_TAG;

		MPI_Recv(&localSum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, DATA_TAG, MPI_COMM_WORLD, &status);
		source = status.MPI_SOURCE;
		offset = getNextTask(task, &allPoints, offset);
		
		if (status.MPI_TAG != TERMINATION_TAG)
		{
			MPI_Send(task, TASK_SIZE, MPI_INT, source, tag, MPI_COMM_WORLD);
		}
		
		globalSum += localSum;
	}

	//Gets the last answer from each proccess:
	for (id = 1; id < numOfProcesses; id++)
	{
		double localSum = 0;
		MPI_Recv(&localSum, 1, MPI_DOUBLE, id, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		globalSum += localSum;
	}

	t1 = MPI_Wtime(); //stop measuring time.

	printf("\nTotal time: %lf\n", t1 - t0);
	printf("\nThe global sum is: %e\n", globalSum);

}

//-------------------Slave---------------------//
void slaveProcess()
{
	int tag;
	int buffer[TASK_SIZE]; 	//buffer will be filled with an X value followed by a Y value
	double mySum = 0;
	MPI_Status status;
	do
	{
		MPI_Recv(&buffer, TASK_SIZE, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		tag = status.MPI_TAG;
		mySum = massive(buffer[0], buffer[1]);

		MPI_Send(&mySum, 1, MPI_DOUBLE, ROOT, tag, MPI_COMM_WORLD);

	} while (tag != TERMINATION_TAG);

}

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
	int i, j;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++) {
			arrX[j + i * size] = i;
			arrY[j + i * size] = j;
		}
	}
}

void checkAllocation(const void* p)
{
	if (!p) {
		printf("Allocation failed.\n");
		exit(0);
	}
	return;
}

//loads into task X and Y values based on offset, and update offset.
int getNextTask(int* task, Points  *allPoints, int offset)
{
	int i;
	for (i = 0; i < TASK_SIZE; i += 2)
	{
		task[i] = allPoints->xArr[offset];
		task[i + 1] = allPoints->yArr[offset];
		offset++;
	}
	return offset;
}