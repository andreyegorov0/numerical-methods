#include <iostream>
#include <random>
#include <iomanip>
#include <fstream>
#include<cstdlib>
#include <cmath>
using namespace std;


double funfunc(double k){   
        return sin(k);
    }

int main(){
    double a = -6.5;
    double b = 6.5;
    double disT;
    disT =  600;
    int N = 30; //количество точек интерполяции
    double* arrX;
    double* arrZ;
    double* arrY;
    double* arrNX;
    double* arrNY;
    srand(time(0)); 
    int c = rand();
    mt19937 randd(c);
    uniform_real_distribution <double> unidistrib(a,b); /// функция создающая равномерное распределние
    arrX = (double*)malloc(N*sizeof(double));
    arrX[0] = a;
    arrX[N-1]= b;

    for (int i = 1; i <= N-2 ; i++ ){ 
         arrX[i] = unidistrib(randd);
        }
    double srav = 0.0;
    for (int i = 1; i < N-1; i ++){
        srav = arrX[i];
        for(int j = 0; j < i; j++){
            if(fabs(srav-arrX[j])<0.08){
                    srav = unidistrib(randd);
                    j=1;
                }
            }
            arrX[i] = srav;
    }



 
    for (int step = 0; step < N -1; ++step) {
        for (int i = 0; i < N - step - 1; ++i) {
            if (arrX[i] > arrX[i + 1]) {
                double temp = arrX[i];
                arrX[i] = arrX[i + 1];
                arrX[i + 1] = temp;
            }
        }
    }
   
    ofstream outfile("data.txt");

    if (!outfile.is_open()){
        cout<< "Error! File is not open!\n";
    }

    for (int j = 0;j <= N-1; j++ ){
        outfile << arrX[j] << " ";

    }

    outfile << "\n";

    arrY = (double*)malloc(N*sizeof(double));

    for (int e = 0; e <= N-1; e++){
        arrY[e] = funfunc(arrX[e]);
    }


    for (int j = 0;j <= N-1;j++ ){
        outfile << arrY[j] << " ";
    }

    arrZ = (double*)malloc((N)*sizeof(double)); /// знаменатель для многочлена Лагранжа
    double itt=1.0;
    for (int i = 0; i<= N-1; i++ ){
        itt=1.0;
        for (int j = 0; j<=N-1; j++ ){
            if (i!=j){
                itt *= (arrX[i] - arrX[j]);
                }
            }
            arrZ[i] = itt; 
        } 

    arrNX = (double*)malloc(disT*sizeof(double));
    arrNY = (double*)malloc(disT*sizeof(double));
    for (int i = 0; i <= disT-1; i++){
        arrNX[i] = a + ((b-a)*i)/(disT-1);
    }
    double s = 0.0;
    double p = 1.0;
    for (int m = 0; m < disT; m++){
            s = 0.0;
            for (int i = 0; i< N; i++){
                p = 1.0;
                for (int j = 0; j< N; j++){
                    if (i!=j){
                      p *= (arrNX[m]-arrX[j]);
                    }
                    
                }
                s += (arrY[i]*p)/(arrZ[i]);


             }
            arrNY[m] = s;


    }
    for (int u = 0; u < N; u++ ){
      cout << arrZ[u] << " ";

    }
    outfile << "\n";
    for (int u = 0; u <= disT-1; u++ ){
        outfile << arrNX[u] << " ";

    }

    outfile << "\n";
    for (int u = 0; u <= disT-1; u++) {

        outfile << arrNY[u] << " ";
    }

   
    outfile.close();

    system ("python3 grafics.py");
    free(arrX);
    free(arrZ);
    free(arrY);
    free(arrNX);
    free(arrNY);
   return 0;

}
 