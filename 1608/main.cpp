#include <iostream>
#include <mpi.h>
#include <math.h>

#define MAX_NZ 10000
#define EPS 1e-8

using namespace std;

int ProcNum = 0;
int ProcRank = 0;
size_t nzA, nzB;

struct crsMatrix {
	double *val;        // Массив значений элементов
	int *col;           // Массив номеров столбцов
	int *rowIndex;		// Массив индексов строк
	size_t nz;          // Количество ненулевых элементов
	size_t size;        // Размер матрицы
};

void parallelMult(crsMatrix *res, crsMatrix A, crsMatrix B);

void initMatrix(crsMatrix *mtx, size_t nz, size_t size) {
	mtx->nz = nz;
	mtx->size = size;
	mtx->val = new double[nz];
	mtx->col = new int[nz];
	mtx->rowIndex = new int[size + 1];
}

void freeMatrix(crsMatrix* mtx) {
	delete[] mtx->val;
	delete[] mtx->col;
	delete[] mtx->rowIndex;
}

void transposition(crsMatrix *AT, crsMatrix A) {
	initMatrix(AT, A.nz, A.size);
	memset(AT->rowIndex, 0, (A.size + 1) * sizeof(int));
	int j1, j2;
	int index;

	for (int i = 0; i < A.nz; i++) {
		AT->rowIndex[A.col[i] + 1]++;
	}

	int sum = 0;

	for (int i = 1; i < A.size + 1; i++) {
		index = AT->rowIndex[i];
		AT->rowIndex[i] = sum;
		sum += index;
	}

	for (int i = 0; i < A.size; i++) {
		j1 = A.rowIndex[i], j2 = A.rowIndex[i + 1];
		for (int j = j1; j < j2; j++) {
			index = AT->rowIndex[A.col[j] + 1];
			AT->val[index] = A.val[j];
			AT->col[index] = i;
			AT->rowIndex[A.col[j] + 1] ++;
		}
	}
}

bool isCompare(crsMatrix A, crsMatrix B) {
	if ((A.nz != B.nz) || (A.size != B.size))
		return false;
	for (int i = 0; i < A.nz; i++)
		if ((A.col[i] != B.col[i]) || (fabs(A.val[i] - B.val[i]) > EPS))
			return false;
	return true;
}

void Multiplicate(crsMatrix *res, crsMatrix A, crsMatrix B) {
	size_t size = A.size;
	size_t nZ = 0;
	int bCol[MAX_NZ];
	double bVal[MAX_NZ];
	int* bRowIndex = new int[size + 1];
	int iBuff = 0, jBuff = 0;
	int col, index;
	double sum;
	int *tmp = new int[size];

	for (int i = 0; i < size; i++) {
		memset(tmp, -1, size * sizeof(int));

		for (int j = B.rowIndex[i]; j < B.rowIndex[i + 1]; j++) {
			col = B.col[j];
			tmp[col] = j;
		}

		for (int j = 0; j < size; j++) {
			sum = 0.0;
			for (int k = A.rowIndex[j]; k < A.rowIndex[j + 1]; k++) {
				col = A.col[k];
				index = tmp[col];
				if (index != -1)
					sum += A.val[k] * B.val[index];
			}
			if (fabs(sum) > EPS) {
				nZ++;
				bCol[iBuff] = j;
				bVal[iBuff] = sum;
				iBuff++;
			}
		}
		bRowIndex[i + 1] = nZ;
	}
	free(tmp);
	free(bRowIndex);
	cout << nZ;

	initMatrix(res, nZ, size);

	for (int i = 0; i < nZ; i++) {
		res->col[i] = bCol[i];
		res->val[i] = bVal[i];
	}

	for (int i = 0; i < size + 1; i++) {
		res->rowIndex[i] = bRowIndex[i];
	}
}

void feelSpMat(crsMatrix *mtx, double *val, int *col, int *rowIndex) {
	for (int i = 0; i < mtx->nz; i++) {
		mtx->val[i] = val[i];
		mtx->col[i] = col[i];
	}
	for (int i = 0; i < mtx->size + 1; i++)
		mtx->rowIndex[i] = rowIndex[i];
}

double indexMatrix(crsMatrix mtx, int i, int j) {
	double res = 0.;
	int N1 = mtx.rowIndex[j];
	int N2 = mtx.rowIndex[j + 1];
	for (int k = N1; k < N2; k++) {
		if (mtx.col[k] == i) {
			res = mtx.val[k];
			break;
		}
	}
	return res;
}

void printMatrix(int size, crsMatrix mtx) {
	//for (int i = 0; i < mtx.size; i++) {
	//	for (int j = 0; j < mtx.size; j++)
	//		cout << indexMatrix(mtx, i, j) << "\t";
	//	cout << endl;
	//}
	int k = mtx.nz;
	cout << "Value: ";
	for (int i = 0; i < k; i++)
		cout << mtx.val[i] << " ";
	cout << endl << "Col:   ";
	for (int i = 0; i < k; i++)
		cout << mtx.col[i] << " ";

	cout << endl << "RowIndex: ";
	for (int i = 0; i < size + 1; i++)
		cout << mtx.rowIndex[i] << " ";
	cout << endl;
}

void generateMatrix(crsMatrix *mtx, int size, int cntInRow) {
	int f, col;
	int notNull = cntInRow * size;
	initMatrix(mtx, notNull, size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < cntInRow; j++) {
			do {
				col = rand() % size;
				f = 0;
				for (int k = 0; k < col; k++)
					if (col == mtx->col[i*cntInRow + k])
						f = 1;
			} while (f == 1);
			mtx->col[i*cntInRow + j] = col;
		}
		for (int j = 0; j < cntInRow - 1; j++)
			for (int k = 0; k < cntInRow - 1; k++)
				if (mtx->col[i*cntInRow + k] > mtx->col[i*cntInRow + k + 1]) {
					int tmp = mtx->col[i*cntInRow + k];
					mtx->col[i*cntInRow + k] = mtx->col[i*cntInRow + k + 1];
					mtx->col[i*cntInRow + k + 1] = tmp;
				}
	}
	for (int i = 0; i < notNull; i++)
		mtx->val[i] = (rand() % 9) / 10.0 + 1;
	int c = 0;
	for (int i = 0; i <= size; i++) {
		mtx->rowIndex[i] = c;
		c += cntInRow;
	}
}

int main(int argc, char **argv) {
	int size = (argc != 1) ? atoi(argv[1]) : 5;
	int cntInRow = (argc != 1) ? atoi(argv[2]) : 2;
	double t1, t2;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	crsMatrix AT, A, B;

	if (ProcRank == 0) {
		generateMatrix(&A, size, cntInRow);
		generateMatrix(&B, size, cntInRow);

		nzA = A.nz;
		nzB = B.nz;
		size = A.size;
		transposition(&AT, A);
	}
	MPI_Bcast(&nzA, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nzB, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

	if (ProcRank != 0) {
		initMatrix(&A, nzA, size);
		initMatrix(&B, nzB, size);
	}

	MPI_Bcast(A.val, nzA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(A.col, nzA, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(A.rowIndex, size + 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(B.val, nzB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(B.col, nzB, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(B.rowIndex, size + 1, MPI_INT, 0, MPI_COMM_WORLD);

	crsMatrix C;
	if (ProcRank == 0) {
		t1 = MPI_Wtime();
	}

	parallelMult(&C, A, B);

	if (ProcRank == 0) {
		t2 = MPI_Wtime();
		cout << "Matrix A:" << endl;
		printMatrix(size, AT);
		cout << endl << "Matrix B:" << endl;
		printMatrix(size, B);
		cout << endl << "Matrix C:" << endl;
		printMatrix(size, C);
		cout << endl << "Time = " << t2 - t1;

		freeMatrix(&C);
		freeMatrix(&AT);
	}
	freeMatrix(&A);
	freeMatrix(&B);

	MPI_Finalize();
	return 0;
}

void parallelMult(crsMatrix *res, crsMatrix A, crsMatrix B) {
	size_t size = A.size;
	size_t NZ = 0;
	int partSize = (int)size / ProcNum;
	int *bCol = new int[MAX_NZ];
	double *bVal = new double[MAX_NZ];
	int *bRowIndex = new int[partSize];
	int iBuff = 0;
	int col, index;
	double sum;
	int *tmp = new int[size];

	bRowIndex[0] = 0;
	for (int i = ProcRank * partSize; i < (ProcRank + 1) * partSize; i++) {
		memset(tmp, -1, size * sizeof(int));

		for (int j = B.rowIndex[i]; j < B.rowIndex[i + 1]; j++) {
			col = B.col[j];
			tmp[col] = j;
		}

		for (int j = 0; j < size; j++) {
			sum = 0.0;
			for (int k = A.rowIndex[j]; k < A.rowIndex[j + 1]; k++) {
				col = A.col[k];
				index = tmp[col];

				if (index != -1) {
					sum += A.val[k] * B.val[index];
				}
			}

			if (fabs(sum) > EPS) {
				NZ++;
				bCol[iBuff] = j;
				bVal[iBuff] = sum;
				iBuff++;
			}

		}
		bRowIndex[i + 1 - ProcRank * partSize] = NZ;
	}

	int *displs = new int[ProcNum];
	int *rcounts = rcounts = new int[ProcNum];

	MPI_Gather(&NZ, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	size_t nZ;

	if (ProcRank == 0) {

		nZ = rcounts[0];
		displs[0] = 0;
		for (int i = 1; i < ProcNum; i++) {
			nZ += rcounts[i];
			displs[i] = displs[i - 1] + rcounts[i - 1];
		}
		initMatrix(res, nZ, size);
	}

	MPI_Gatherv(bVal, NZ, MPI_DOUBLE, res->val, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(bCol, NZ, MPI_INT, res->col, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(bRowIndex + 1, partSize, MPI_INT, res->rowIndex + 1, partSize, MPI_INT, 0, MPI_COMM_WORLD);

	if (ProcRank == 0) {
		int endIndex = res->rowIndex[partSize];
		for (int i = 1; i < ProcNum; i++) {
			for (int j = 1 + i * partSize; j <= partSize * (i + 1); j++) {
				res->rowIndex[j] += endIndex;
			}
			endIndex = res->rowIndex[partSize * (i + 1)];
		}
		res->rowIndex[0] = 0;
		free(tmp);
	}
}