#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <time.h>

using namespace std;

int main(int argc, char *argv[]) {
    int count; // Количество символов

    cout << "Enter count ";
    cin >> count;

    char* str = NULL;
    int ProcNum; // Число процессов
    int ProcRank; // Номер ранга
    int num = 0;
    double sum = 0, t;
    setlocale(LC_ALL, "Russian");

    MPI_Status stat; // Статус выполнения операции передачи данных
    MPI_Init(&argc, &argv); // Инициализации среды выполнения программы
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum); // Получение количества процессов 
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank); // Определение ранга процесса

    if (ProcRank == 0) { // Заполняем строку и передаем другим процессам

        char *str = new char[count];
        char arr[33];
        srand(time(0));
        for (int i = 0, j = 192; i < 32; i++) {
            arr[i] = (char)j;
            j++;
        }
        arr[32] = ' ';
        for (int i = 0; i < count; i++)
            str[i] = arr[rand() % 33];

        if (ProcNum < 1) {
            printf("Size of string < 1\n");
            return 0;
        }
        else
            cout << "Number proc: " << ProcNum << endl;
    }

    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);       // Передача данных
    MPI_Bcast(&str, count, MPI_CHAR, 0, MPI_COMM_WORLD);

    int size = count / ProcNum;
    int i1 = size * ProcRank;
    int i2 = size * (ProcRank + 1);

    cout << ProcRank << " proc start work\n";
    cout << str << endl;

    if (ProcRank != 0) {
        for (int i = ProcRank; i < count; i += ProcNum) {
            MPI_Recv(str, 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &stat);
            for (int j = i1; j< i2; j++) {
                if (str[i] == ' ')
                    sum++;
                num++;
            }
            MPI_Send(&sum, 1, MPI_INT, 0, i, MPI_COMM_WORLD);

        }
    }
    
    if (ProcRank == 0) {
        for (int i = 0; i < count; i++) {
            if (size != 0)
                MPI_Recv(&sum + i, 1, MPI_DOUBLE, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &stat);
        }
    }

    MPI_Finalize();
    return 0;
}