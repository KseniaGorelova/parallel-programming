#include<mpi.h>
#include <time.h>
#include<iostream>

using namespace std;

void my_sum(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) {
	if ((type == MPI_INT) || (type == MPI_DOUBLE) || (type == MPI_FLOAT)) {
		for (int i = 0; i < count; i++) {
			if (type == MPI_INT)
				((int*)recvbuf)[i] += ((int*)sendbuf)[i];
			if (type == MPI_DOUBLE)
				((double*)recvbuf)[i] += ((double*)sendbuf)[i];
			if (type == MPI_FLOAT)
				((float*)recvbuf)[i] += ((float*)sendbuf)[i];
		}
	}
}

void my_sum_tree(void* sendbuf, void* recvbuf, int count){
	int ProcNum;
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcNum);
	for (int i = 0; i < count; i++) 
		((int*)recvbuf)[i] += ((int*)sendbuf)[i];


}

void my_prod(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) {
	if ((type == MPI_INT) || (type == MPI_DOUBLE) || (type == MPI_FLOAT)) {
		for (int i = 0; i < count; i++) {
			if (type == MPI_INT)
				((int*)recvbuf)[i] *= ((int*)sendbuf)[i];
			if (type == MPI_DOUBLE)
				((double*)recvbuf)[i] *= ((double*)sendbuf)[i];
			if (type == MPI_FLOAT)
				((float*)recvbuf)[i] *= ((float*)sendbuf)[i];
		}
	}
}

void my_min(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) {
	if ((type == MPI_INT) || (type == MPI_DOUBLE) || (type == MPI_FLOAT)) {
		for (int i = 0; i < count; i++) {
			if (type == MPI_INT)
				if (((int*)recvbuf)[i] > ((int*)sendbuf)[i])
					((int*)recvbuf)[i] = ((int*)sendbuf)[i];
			if (type == MPI_DOUBLE)
				if (((double*)recvbuf)[i] > ((double*)sendbuf)[i])
					((double*)recvbuf)[i] = ((double*)sendbuf)[i];
			if (type == MPI_FLOAT)
				if (((float*)recvbuf)[i] > ((float*)sendbuf)[i])
					((float*)recvbuf)[i] = ((float*)sendbuf)[i];
		}
	}
}

void my_max(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) {
	if ((type == MPI_INT) || (type == MPI_DOUBLE) || (type == MPI_FLOAT)) {
		for (int i = 0; i < count; i++) {
			if (type == MPI_INT)
				if (((int*)recvbuf)[i] < ((int*)sendbuf)[i])
					((int*)recvbuf)[i] = ((int*)sendbuf)[i];
			if (type == MPI_DOUBLE)
				if (((double*)recvbuf)[i] < ((double*)sendbuf)[i])
					((double*)recvbuf)[i] = ((double*)sendbuf)[i];
			if (type == MPI_FLOAT)
				if (((float*)recvbuf)[i] < ((float*)sendbuf)[i])
					((float*)recvbuf)[i] = ((float*)sendbuf)[i];
		}
	}
}

int my_reduce(void *sendbuf, void* recvbuf, int count, MPI_Datatype type,
	MPI_Op op, int root, MPI_Comm comm) {

	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status status;

	if (ProcRank == root) {
		for (int i = 0; i < count; i++) {
			if (type == MPI_INT)
				((int*)recvbuf)[i] = ((int*)sendbuf)[i];
			if (type == MPI_DOUBLE)
				((double*)recvbuf)[i] = ((double*)sendbuf)[i];
			if (type == MPI_FLOAT)
				((float*)recvbuf)[i] = ((float*)sendbuf)[i];
		}
		MPI_Recv(sendbuf, count, type, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
		if (op == MPI_SUM)
			my_sum(sendbuf, recvbuf, count, type);
		if (op == MPI_PROD)
			my_prod(sendbuf, recvbuf, count, type);
		if (op == MPI_MIN)
			my_min(sendbuf, recvbuf, count, type);
		if (op == MPI_MAX)
			my_max(sendbuf, recvbuf, count, type);
	}
	else
		MPI_Send(sendbuf, count, type, root, 0, comm);
}

void tree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type,
	MPI_Op op, int root, MPI_Comm comm, int* sendmas, int size, int h) {
	int curSize, ProcNum, ProcRank;
	int* curMas;
	int* elem = new int[count];

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status status;

	if (size % 2) {
		curSize = size / 2 + 1;
		curMas = new int[curSize];
	}
	else {
		curSize = size / 2;
		curMas = new int[curSize];
	}
	for (int i = 0; i < curSize; i++)
		curMas[i] = -1;
	h++;
	for (int i = 0, j = 0; i < size; i++) {
		if (sendmas[i] == ProcRank) {
			if (!(size % 2)) {
				if (!(i % 2)) {
					MPI_Send(sendbuf, count, type, sendmas[i + 1], 0, comm);
					for (int k = 0; k < count; k++)
						elem[k] = ((int *)(recvbuf))[k];
				}
				else {
					MPI_Recv(recvbuf, count, type, sendmas[i - 1], 0, comm, &status;
					my_sum_tree(sendbuf, recvbuf, count);
					sendbuf = recvbuf;
				}
			}
			else {
				if ((ProcRank != ProcNum - 1)) {
					if (!(i % 2)) {
						MPI_Send(sendbuf, count, type, sendmas[i + 1], 0, comm);
						for (int k = 0; k < count; k++)
							elem[k] = ((int*)(recvbuf))[k];
					}
					else {
						for (int k = 0; k < count; k++) {
							((int*)recvbuf)[i] = ((int*)sendbuf)[i];
						}
						MPI_Recv(sendbuf, count, type, sendmas[i - 1], 0, comm, &status);
						my_sum_tree(sendbuf, recvbuf, count);
						for (int k = 0; k < count; k++) {
							((int*)sendbuf)[i] = ((int*)recvbuf)[i];
						}
					}
				}
				else {
					curMas[j] = i;
					j++;
				}
			}
		}
		if (!(size % 2)) {
			if (i % 2) {
				curMas[j] = i;
				j++;
			}
		}
	}
	if (curSize != 1)
		tree(sendbuf, elem, count, type, op, root, comm, curMas, curSize, h);
	else {
		if (curMas[0] != root) {
			if (ProcRank == curMas[0])
				MPI_Send(sendbuf, count, type, root, 0, comm);
			if (ProcRank == root)
				MPI_Recv(recvbuf, count, type, curMas[0], 0, comm, &status);
		}
		else
			recvbuf = sendbuf;
	}
}

int tree_reduce(void *sendbuf, void *recvbuf, int count,
	MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	int ProcNum, ProcRank;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int* sendmas = new int[ProcNum];
	int recvmas;
	for (int i = 0; i < ProcNum; i++)
		sendmas[i] = i;
	tree(sendbuf, recvbuf, count, type, op, root, comm, sendmas, ProcNum, 0);
	return 0;
}

int main(int argc, char *argv[]) {
	int ProcRank, ProcNum;

	srand(1);
	double t1, tTree1, tMPI1, t2, tTree2, tMPI2;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int* arr1 = new int[20];
	int* arr2 = new int[20];

	for (int i = 0; i < 20; i++)
		arr1[i] = rand() % 100;

	if (ProcRank == 0)
		t1 = MPI_Wtime();

	my_reduce(arr1, arr2, 20, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0) {
		t2 = MPI_Wtime();
		tTree1 = MPI_Wtime();
	}
	tree_reduce(arr1, arr2, 20, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0) {
		tTree2 = MPI_Wtime();
		tMPI1 = MPI_Wtime();
	}
	MPI_Reduce(arr1, arr2, 20, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0) {
		tMPI2 = MPI_Wtime();
		cout << "Time = " << t2 - t1 << endl;
		cout << "Time tree = " << tTree2 - tTree1 << endl;
		cout << "Time MPI = " << tMPI2 - tMPI1 << endl;
	}

	MPI_Finalize();
	return 0;
}
