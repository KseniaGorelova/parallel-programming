#define _CRT_SECURE_NO_WARNINGS

#include <mpi.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include<string.h>
#include <iostream>

using namespace std;
//Сделать рандом

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

	//Вводить количество символов с клавиатуры

    int count; // Количество символов
	cout << "Enter count ";
	cin >> count;
    char* str = NULL;
    int ProcNum; // Число процессов
    int ProcRank; // Номер ранга
    int size, i = 0, i1, i2, num = 0;
    double sum = 0;

    MPI_Status stat; // Статус выполнения операции передачи данных
    MPI_Init(&argc, &argv); // Инициализации среды выполнения программы
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum); // Получение количества процессов 
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank); // Определение ранга процесса

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

    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);       // Передача данных
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
			//Переписать на reduce
            MPI_Recv(&sum, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &stat); // Прием данных
            num += sum;
        }


        cout<< str<<endl;
        cout<<"Number of words: "<< num;
    }
    else
        MPI_Send(&sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD); // Передача данных
    MPI_Finalize();
    return 0;
}