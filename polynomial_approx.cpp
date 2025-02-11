#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <random>
#include <cmath>
#include <ctime>
#include </Users/andrejegorov/Desktop/source_dir/Eigen/Sparse>
#include </Users/andrejegorov/Desktop/source_dir/Eigen/SparseCore>
using namespace std;

double Func(double x){ // Обьявление функции
    return cos(x)+sin(x);
}

double phi(double x, int i, int N, double* arrX){ //считает базисные функции в точке х
    double t = 1;
    for(int j=0;j<N;j++) {
        if (j!=i){
            t*=(x-arrX[j])/(arrX[i]-arrX[j]);
        }
    }
    return t;
}

double mul(double* x, int len, double* y){ //скалярное произведение векторов
    double t = 0;
     for (int i=0;i<len;i++){
        t+=x[i]*y[i];
    }
    return t;
}

void selectionSort(double arr[], int n) { //сортировка
    for (int i = 0; i < n - 1; i++) {
        int min = i;
        for (int j = i + 1; j < n; j++) {
            if (arr[j] < arr[min])
                min = j;
        }
       if (min != i) {
            double temp = arr[min];
            arr[min] = arr[i];
            arr[i] = temp;
        }
    }
}


int main(){
    double a = -6.28;
    double b = 6.28;
    int K = 10; // количество конечных элементов(интервалов)
    int N = 4; // количество узлов на каждом конечном элементе
    int M =  K*N-(K-1);//общее количество узлов на отрезке
    int L = 6; //количество случайных точек на интервале
    int N1 = 101; //количество узлов конечного элемента на Мвиз
    double h = (b-a)/(M-1); //h - шаг сетки
    double h1; // шаг мелкой равномерной сетки
    double srav=0;



    double* arrH; // массив нашей равномерной сетки
    double** arrR; //двумерный массив случайных точек
    double* arrOtr; //концы отрезков
    double*** arrPHI_k; //трехмерный массив базисных функий, первый индекс - номер конечного элемента, второй - номер случайной точки на интервале, третий - номер базисной функции на интервале
    double* temp; // просто временный массив
    double** arrPhi; // матрица фи
    double** arrPhi_T; // матрица фи транспонированная
    double** arrFTF; //матричное произведение фи_транспонированная*фи
    double* arrF; // вектор f в случайных точках
    double** arrF2d; //тоже вектор f в случайных точках, но двумерный - до сложения пересечений
    double* arrPhi_F; // фи*f - правая часть слау
    double** arrP; // вектор значений энп
    double*arrFH; // вектор значений функции на равномерной сетке
    double** C_2d; // вектор констант решения С, преобразованный в двумерный массив
    double*** arrPhi_k; // матрица значений базисных функций на мелкой сетке 
    double** arrH1; // матрица мелкой равномерной сетки


    //////////////выделение памяти///////////////////
    
    temp = (double*)malloc(N*sizeof(double));
    arrH = (double*)malloc(M*sizeof(double));
    arrP = (double**)malloc(K*sizeof(double));
    arrFH = (double*)malloc(M*sizeof(double));
    arrOtr = (double*)malloc((K+1)*sizeof(double));
    arrF = (double*)malloc((L*K-K+1)*sizeof(double));
    arrPhi_F = (double*)malloc((N*K-K+1)*sizeof(double));
    arrR = (double**)malloc(K*sizeof(double*));
    arrF2d = (double**)malloc(K*sizeof(double*));
    arrPhi = (double**)malloc((L*K-K+1)*sizeof(double*));
    arrPhi_T = (double**)malloc((N*K-K+1)*sizeof(double*));
    arrFTF = (double**)malloc((N*K-K+1)*sizeof(double*));
    arrPHI_k = (double***)malloc(K*sizeof(double**));
    C_2d = (double**)malloc(K*sizeof(double*));
    arrPhi_k = (double***)malloc(K*sizeof(double**));
    arrH1 = (double**)malloc(K*sizeof(double*));

    

    for (int i = 0;i<K;i++){
        arrR[i] = (double*)malloc(L*sizeof(double));//двумерный массив случайных точек на отрезках
    }

    for (int i=0;i<K;i++){
        arrP[i] = (double*)malloc(N1*sizeof(double));
    }

    for (int i=0;i<K;i++){
        C_2d[i] = (double*)malloc(N*sizeof(double));
    } 

    for (int i=0;i<K;i++){
        arrH1[i] = (double*)malloc(N1*sizeof(double));
    } 


    for (int i = 0;i<K;i++){
        arrF2d[i] = (double*)malloc(L*sizeof(double));
    }

    for (int i=0;i<L*K-K+1;i++){
        arrPhi[i] = (double*)malloc((N*K-K+1)*sizeof(double));
    } 

    for (int i=0;i<(N*K-K+1);i++){
        arrPhi_T[i] = (double*)malloc((L*K-K+1)*sizeof(double));
    } 

    for (int i=0;i<N*K-K+1;i++){
        arrFTF[i] = (double*)malloc((N*K-K+1)*sizeof(double));
    } 

    for (int i=0;i<K;i++){
        arrPHI_k[i] = (double**)malloc(L*sizeof(double*));
        for (int j=0;j<L;j++){
            arrPHI_k[i][j] = (double*)malloc(N*sizeof(double));
        }
    }

    for (int i=0;i<K;i++){
        arrPhi_k[i] = (double**)malloc(N1*sizeof(double*));
        for (int j=0;j<N1;j++){
            arrPhi_k[i][j] = (double*)malloc(N*sizeof(double));
        }
    }

    ///////////////////////////

    for(int i=0;i<M;i++){   //равномерная сетка на (a,b)
        arrH[i] = a+i*h;
    }
    // for (int i=0;i<M;i++){ //вывод массива сетки
    //     cout<<arrH[i]<<' ';
    // }
    // cout<<'\n';

    for(int i=0;i<K+1;i++){   //концы конечных элементов
        arrOtr[i] = arrH[i*(N-1)];
    }

    // for (int i=0;i<K+1;i++){ //вывод массива концов отрезков
    //     cout<<arrOtr[i]<<' ';
    // }
    // cout<<'\n';


    h1 = (arrOtr[1]-arrOtr[0])/(N1-1); // построение мелкой равномерной сетки
    for (int k=0;k<K;k++){
        for (int i = 0; i<N1;i++){
            arrH1[k][i] = arrOtr[k] + i * h1;
        }
    }

    // for (int k=0;k<K;k++){ //вывод мелкой равномерной сетки
    //     for (int i = 0; i<N1;i++){
    //         cout<<arrH1[k][i]<<' ';
    //     }
    //     cout<<'\n';
    // }


    for (int i=0;i<M;i++){     //значения функции на равномерной сетке
        arrFH[i] = Func(arrH[i]);
    }
    // for (int i=0;i<M;i++){ //вывод значений функции на равномерной сетке
    //     cout<<arrFH[i]<<' ';
    // }
    // cout<<'\n';



    for (int k=0;k<K;k++){ //генерация случайных точек на интервалах
    
        std::random_device genSource;

        std::uniform_real_distribution<> generator(arrOtr[k],arrOtr[k+1]);

        for (int index = 0; index < L; index++)
        {
            arrR[k][index] = generator(genSource);
        }
    }


    

    for (int k=0;k<K;k++){ //сортировка случайных точек 
        selectionSort(arrR[k],L);
    }

    // for (int i=0;i<K;i++){ //вывод двумерного массива случайных точек на конечных элементах
    //     for (int j=0;j<L;j++){
    //         cout<<arrR[i][j]<<' ';
    //     }
    //     cout<<'\n';
    // }

    for (int k=0;k<K;k++){ //заполнение матриц фи_к 
        for (int m=0;m<N;m++){
                    temp[m] = arrH[k*(N-1)+m];
            }
        for (int j=0;j<L;j++){
            for (int i=0;i<N;i++){
                arrPHI_k[k][j][i] = phi(arrR[k][j],i,N,temp);
            }
        }
    }


    // for (int j=0;j<L;j++){ //вывод второго блока - фи_2 
    //     for (int i=0;i<N;i++){
    //         cout<<arrPHI_k[1][j][i]<<' ';
    //     }
    //     cout<<'\n';
    // }

    ///////////////составление матрицы фи из блоков фи_к/////////////////

    for (int k=0;k<K-1;k++){
        for (int j=0;j<L-1;j++){
            for (int i=0;i<N;i++){
                 arrPhi[k*(L-1)+j][k*(N-1)+i]= arrPHI_k[k][j][i];
            }
        }
        for (int i=0;i<2*N-1;i++){
            if (i<N-1){
                arrPhi[(k+1)*(L-1)][k*(N-1)+i]= arrPHI_k[k][(L-1)][i];
            }
            else {
                arrPhi[(k+1)*(L-1)][k*(N-1)+i]= arrPHI_k[k+1][(L-1)][i-N+1];
            }
        }
    }
    for (int i = 0;i<L;i++){
        for (int j=0;j<N;j++){
            arrPhi[(K-1)*(L-1)+i][(K-1)*(N-1)+j]= arrPHI_k[K-1][i][j];
        }
        
    } 
    for (int k=0;k<K-1;k++){
        arrPhi[(k+1)*(L-1)][k*(N-1)+N-1]= arrPHI_k[k][(L-1)][N-1]+arrPHI_k[k+1][0][0];
    }


    // ofstream outfile("data.txt"); //запись матрицы фи в файл
    // if (!outfile.is_open()){
    //     cout<<"Error! File is not open!\n";
    // }
    // for (int i=0;i<L*K-K+1;i++){
    //     for (int j=0;j<N*K-K+1;j++){
    //         outfile<<arrPhi[i][j]<<'\t';
    //     }
    //     outfile<<'\n';
    // }

    // outfile.close();


    /////////////////////////////////////


    for (int i=0;i<L*K-K+1;i++){ //матрица фи_транспонированная
        for (int j=0;j<N*K-K+1;j++){
            arrPhi_T[j][i]=arrPhi[i][j];
        }
    }

    // ofstream outfile1("data1.txt"); //запись матрицы фи_т в файл
    // if (!outfile1.is_open()){
    //     cout<<"Error! File is not open!\n";
    // }
    // for (int i=0;i<N*K-K+1;i++){
    //     for (int j=0;j<L*K-K+1;j++){
    //         outfile1<<arrPhi_T[i][j]<<'\t';
    //     }
    //     outfile1<<'\n';
    // }
    // outfile1.close();


    for (int i=0;i<N*K-K+1;i++){ //матричное произведение фи_транспонированная*фи
        for (int j=0;j<N*K-K+1;j++){
            arrFTF[i][j]=mul(arrPhi_T[i],L*K-K+1,arrPhi_T[j]);
        }
    }

    // ofstream outfile2("data2.txt"); //запись матрицы фи_т*фи в файл
    // if (!outfile2.is_open()){
    //     cout<<"Error! File is not open!\n";
    // }
    // for (int i=0;i<N*K-K+1;i++){
    //     for (int j=0;j<N*K-K+1;j++){
    //         outfile2<<arrFTF[i][j]<<'\t';
    //     }
    //     outfile2<<'\n';
    // }
    // outfile2.close();


    for (int k=0;k<K;k++){  //заполнение двумерного вектора f в случайных точках
        for (int i=0;i<L;i++){
            arrF2d[k][i] = Func(arrR[k][i]);
        }
    }

    // for (int i=0;i<K;i++){ //вывод двумерного вектора значений f в случайных точках 
    //     for (int j=0;j<L;j++){
    //         cout<<arrF2d[i][j]<<'\t';
    //     }
    //     cout<<'\n';
    // }

    // ofstream outfile3("data3.txt"); //запись двумерного вектора значений f в случайных точках в файл
    // if (!outfile3.is_open()){
    //     cout<<"Error! File is not open!\n";
    // }
    // for (int i=0;i<K;i++){
    //     for (int j=0;j<L;j++){
    //         outfile3<<arrF2d[i][j]<<'\t';
    //     }
    //     outfile3<<'\n';
    // }

    // outfile3.close();

    for (int k=0;k<K-1;k++){ 
        arrF2d[k][L-1]+=arrF2d[k+1][0];
    }

    
    // for (int i=0;i<K;i++){ //снова вывод двумерного вектора значений f в случайных точках, сложили пересекающиеся эл-ты
    //     for (int j=0;j<L;j++){
    //         cout<<arrF2d[i][j]<<'\t';
    //     }
    //     cout<<'\n';
    // }
    
    //преобразование двумерного вектора значений f в одномерный, пересекающиеся элементы складываются//

    arrF[0]=arrF2d[0][0];
    for (int k=0;k<K;k++){
        for (int i=1;i<L;i++){
            arrF[k*(L-1)+i]=arrF2d[k][i];
        }
    }


    // for (int i=0;i<L*K-K+1;i++){ //вывод результирующего вектора значений функции в случайных точках после сложения
    //     cout<<arrF[i]<<' ';
    // }
    // cout<<'\n';

    ////////////////////////


    for (int i=0;i<N*K-K+1;i++){ //вычисление правой части матричного уравнения - фи_т*f
        arrPhi_F[i] = mul(arrPhi_T[i],(L*K-K+1),arrF);
    }

    //  for (int i=0;i<N*K-K+1;i++){ //вывод вектора значений правой части
    //     cout<<arrPhi_F[i]<<' ';
    // }
    // cout<<'\n';


    ////////////решение системы библиотекой питона для проверки результатов///////////////

    ofstream out_A("A.txt");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            out_A << arrFTF[i][j] << " ";
        }
        out_A << endl;
    }
    out_A.close();
    //Запись вектора в файл
    ofstream out_B("B.txt");
    for (int i = 0; i < M; i++) {
        out_B << arrPhi_F[i] << endl;
    }
    out_B.close();
    //Решаем СЛАУ внутри питона
    system("python3 solve.py");

    ifstream input_coefs("X.txt");  //файл для чтения
    double* C = new double[M];  //вектор-решение СЛАУ - наши коэффициенты при базисных функциях
    int i = 0;
    double c;
    //считываем файл с решением
    while (input_coefs >> c && i < M) {
        C[i] = c;
        i++;
    }


    

    //////////////////////////////

    ////////////решение слау библиотекой eigen//////////////////
    Eigen::SparseMatrix<double> A(M, M);
    Eigen::VectorXd B(M);
    Eigen::VectorXd X(M);
    for (int i = 0; i < M; i++) { B[i] = arrPhi_F[i]; }
    std::vector<Eigen::Triplet<double> > coefficients;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            coefficients.push_back(Eigen::Triplet<double>(i, j, arrFTF[i][j]));
        }
    }

    A.setFromTriplets(coefficients.begin(), coefficients.end());
    Eigen::ConjugateGradient<Eigen::SparseMatrix < double >,Eigen::Lower|Eigen::Upper> cg;
    cg.compute(A);
    X = cg.solve(B);
    // for (int i=0;i<M;i++){
    //     cout<<X[i]<<' '; 
    // }
    for (int i=0;i<M;i++){
        C[i]=X[i]; //в векторе С получили искомый набор коэффициентов
    }

    /////////////////////////////////////

    ////////преобразование С в двумерный вектор для получения энп///////

    // for (int i=0;i<M;i++){ // вывод вектора С до преобразования
    //     cout<<C[i]<<' ';
    // }
    // cout<<'\n'<<'\n';

    C_2d[0][0]=C[0]; //преобразование вектора решения С в двумерный массив
    for (int k=0;k<K;k++){
        for (int i=1;i<N;i++){
            C_2d[k][i]=C[k*(N-1)+i];
        }
        if (k!=0){
            C_2d[k][0]=C[k*(N-1)];
        }
    }
  
    // for (int k=0;k<K;k++){ // вывод двумерного преобразованного вектора С
    //     for (int i=0;i<N;i++){
    //         cout<<C_2d[k][i]<<" ";
    //     }
    //     cout<<'\n';
    // }



    for (int k=0;k<K;k++){ //заполнение матриц фи_к новой для мелкой равномерной сетки (вектор С остается старый, посчитанный по случайным точкам)
        for (int m=0;m<N;m++){
                    temp[m] = arrH[k*(N-1)+m];
            }
        for (int j=0;j<N1;j++){
            for (int i=0;i<N;i++){
                arrPhi_k[k][j][i] = phi(arrH1[k][j],i,N,temp);
            }
        }
    }

    for (int k=0;k<K;k++){ // построение энп по старому вектору С и новой матрице фи(она пересчитывается для мелкой равномерной сетки)
        for (int j=0;j<N1;j++){
            arrP[k][j]=mul(arrPhi_k[k][j],N,C_2d[k]);
        }
    }



    // for (int i=0;i<M;i++){ //подставляем найденные коэффициенты С - получаем искомый ЭНП для f 
    //     arrP[i] = mul(arrPhi[i],M,C);
    // }

    ////////////вывод результатов в файл для графики питона////////////////

    ofstream out_p("p.txt");
    for (int k = 0; k < K; k++) {
        for (int j=0;j<N1;j++){
            out_p << arrP[k][j]<<' ';
        }
        out_p<<'\n';
    }
    out_p.close();

    ofstream out_h("h.txt");
    for (int k = 0; k < K; k++) {
        for (int j=0;j<N1;j++){
            out_h <<arrH1[k][j]<<' ';
        }
        out_h<<'\n';
    }
    out_h.close();

    //////////////////////////////////


    ////////////промежуточные выводы для отладки/////////////


    // for (int i=0;i<K;i++){
    //     for (int j=0;j<L;j++){
    //         cout<<arrR[i][j]<<endl;
    //     }
    //     cout<<'\n';
    // }

    //  for (int j=0;j<L;j++){
    //         for (int i=0;i<N;i++){
    //             cout<<arrPHI_k[0][j][i]<<' ';
    //         }
    //         cout<<'\n';
    //     }
    // cout<<'\n';

    // for (int j=0;j<L;j++){
    //         for (int i=0;i<N;i++){
    //             cout<<arrPHI_k[1][j][i]<<' ';
    //         }
    //         cout<<'\n';
    //     }
    // cout<<'\n';

    // for (int j=0;j<L;j++){
    //         for (int i=0;i<N;i++){
    //             cout<<arrPHI_k[2][j][i]<<' ';
    //         }
    //         cout<<'\n';
    //     }
    // cout<<'\n';

    // for (int j=0;j<L;j++){
    //         for (int i=0;i<N;i++){
    //             cout<<arrPHI_k[8][j][i]<<' ';
    //         }
    //         cout<<'\n';
    //     }
    // cout<<'\n';

    // for (int j=0;j<L;j++){
    //         for (int i=0;i<N;i++){
    //             cout<<arrPHI_k[9][j][i]<<' ';
    //         }
    //         cout<<'\n';
    //     }

    // for (int j=0;j<L;j++){
    //         for (int i=0;i<N;i++){
    //             cout<<arrPHI_k[10][j][i]<<' ';
    //         }
    //         cout<<'\n';
    //     }

    // for (int k=0;k<K;k++){
    //     for (int m=0;m<N;m++){
    //         temp[m] = arrH[k*(N-1)+m];
    //     }
    //     for (int b=0;b<N;b++){
    //         cout << temp[b] << " ";
    //     }
    //         cout<<'\n'<<endl;
    // }


    //////////////////////////////////


    system ("python3 graf.py");//графика в питоне



    //////////////очистка памяти///////////////////

    

    
    for (int i=0;i<K;i++){
        for (int j=0;j<L;j++){
            free(arrPHI_k[i][j]);
        }
    }
    for (int i=0;i<K;i++){
       free(arrPHI_k[i]);
    }

    for (int i=0;i<K;i++){
        for (int j=0;j<N1;j++){
            free(arrPhi_k[i][j]);
        }
    }
    for (int i=0;i<K;i++){
       free(arrPhi_k[i]);
    }

    for (int i = 0;i<K;i++){
        free(arrR[i]);
    }

    for (int i=0;i<L*K-K+1;i++){
        free(arrPhi[i]);
    }
    for (int i=0;i<N*K-K+1;i++){
        free(arrPhi_T[i]);
    }
    for (int i=0;i<N*K-K+1;i++){
        free(arrFTF[i]);
    }
    for (int k=0;k<K;k++){
        free(arrF2d[k]);
    }
    for (int k=0;k<K;k++){
        free(C_2d[k]);
    }
    for (int k=0;k<K;k++){
        free(arrH1[k]);
    }

    free(C_2d);
    free(arrF2d);
    free(arrPhi); 
    free(arrPhi_T);
    free(arrFTF);
    free(arrPHI_k);
    free(arrR);
    free(arrH);
    free(arrOtr);
    free(temp);
    free(arrF);
    free(arrFH);
    free(arrP);
    free(arrPhi_F);
    free(arrPhi_k);
    free(arrH1);
    delete[] C;

    ////////////////////////////


    return 0;
}