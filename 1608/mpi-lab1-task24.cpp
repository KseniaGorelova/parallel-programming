#define _CRT_SECURE_NO_WARNINGS

#include <mpi.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include<string.h>
#include <iostream>

using namespace std;
//������� ������

//char * FillStr(int size) {
//    
//    //char arr[33];
//    //int i = 0, j;
//    //for (i = 0, j = 192; i < 32; i++) {
//    //    arr[i] = (char)j;
//    //    j++;
//    //}
//    //arr[32] = ' ';
//    
//    return str;
//}

int main(int argc, char *argv[]) {

	//������� ���������� �������� � ����������

    int count; // ���������� ��������
	cout << "Enter count ";
	cin >> count;
    char* str = NULL;
    int ProcNum; // ����� ���������
    int ProcRank; // ����� �����
    int size, i = 0, i1, i2, num = 0;
    double sum = 0;

    MPI_Status stat; // ������ ���������� �������� �������� ������
    MPI_Init(&argc, &argv); // ������������� ����� ���������� ���������
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum); // ��������� ���������� ��������� 
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank); // ����������� ����� ��������

    if (ProcRank == 0) {
		str = (char*)malloc(count);
		cin >> str;

		if (ProcNum < 1) {
			cout<<"Size of string < 1\n";
			return 0;
		}
		else
			cout << "Number proc: " << ProcNum << endl;
    }

    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);       // �������� ������
    MPI_Bcast(&str, count, MPI_CHAR, 0, MPI_COMM_WORLD);

    size = count / ProcNum;
    i1 = size * ProcRank;
    i2 = size * (ProcRank + 1);

    cout<< ProcRank<< " proc start work\n";

    for (i = i1; i < i2; i++) {
        if (str[i] == ' ')
            sum++;
    }
    num++;

    if (ProcRank == 0) {
        num = sum;
        for (i = 1; i < ProcNum; i++) {
			//���������� �� reduce
            MPI_Recv(&sum, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &stat); // ����� ������
            num += sum;
        }


        cout<< str<<endl;
        cout<<"Number of words: "<< num;
    }
    else
        MPI_Send(&sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD); // �������� ������
    MPI_Finalize();
    return 0;
}