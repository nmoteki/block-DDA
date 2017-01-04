#include<iostream>
#include<omp.h>
using namespace std;

int main(){
    int hh=0;
    #pragma omp parallel for
    for(int i=1; i<=1000; ++i){
        hh+=i;
    }
    cout << "hh= " << hh << endl;

    cout << "omp_get_max_threads()= " << omp_get_max_threads() << endl;
    omp_set_num_threads(4);



int hh=0;
int a[100];
int b[100];
//int i;
#pragma omp parallel for reduction(+:hh)
for(int i=0; i<100; ++i){
   //#pragma omp atomic
   hh+=i;
   a[i]=i;
   b[i]=omp_get_thread_num();
}
cout << "hh= " << hh << endl;

for(int i=0; i<100; ++i){
   cout << "a= " << a[i] << endl;
   cout << "b= " << b[i] << endl;
}



int suma=0;
#pragma omp parallel for reduction(+:suma)
for(int i=1; i<=100; ++i){
   suma+=a[i];
   b[i]= omp_get_thread_num();
}
cout << "suma= " << suma << endl;
for(int i=0; i<100; ++i){
   cout << "b= " << b[i] << endl;
}



/*
int hh=0;
#pragma omp parallel for reduction(+:hh)
for(i=0; i<100; ++i){
   //x=i;
   hh+=i;
}
cout << "hh= " << hh << endl;
*/


/*
    int aa = 0;
#pragma omp parallel reduction(+:aa)
{
   //cout << "omp_get_num_threads() " << omp_get_num_threads() << endl;
   aa++;
}
cout << "aa= " << aa << endl;
*/
}
