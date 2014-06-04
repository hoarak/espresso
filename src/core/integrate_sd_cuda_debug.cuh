#ifndef __helper_functions_cuh__
#define __helper_functions_cuh__

#include "config.hpp"

#include <sys/time.h>
#include <assert.h>
#include <iostream>
#include <iomanip> 

#ifdef SD_USE_FLOAT
typedef float real;
#else
typedef double real;
#endif


// template implementation below
template<typename T> void printMatrixDev( T * data, int m, int n, const char * msg);
template<typename T> void printMatrixDev( T * data, int m, int n);
template<typename T> void printMatrixHost( T * data, int m, int n, const char * msg);

template<typename T> void printDistances( T * r_h, T, int N, T min, T max);
template<typename T> void printDistances( T * r_h, T, int N);

void printPosDev( real * data, int m, const char * msg);
void printVectorDev( real * data, int m);
void printVectorDev( real * data, int m, const char * msg);
//void printVectorDev( float * data, int m);
void printVectorHost( real * data, int m);
void printVectorHost( real * data, int m, const char * msg);
void cudaCheckError(const char *msg);
bool hasAnyNanDev(const real * data, int length);
bool isSymmetricDev(const real * data, int lda, int size);


real getAbsMax(real * data, int length);
real getAbsMaxDev(real * data_d, int length);
int countMinusDev(real * data_d, int length);

class myTimer {
public:
  // init numCounters counter
  myTimer(unsigned int _numCounters);
  ~myTimer();
  // add the time since the last call to counter _counter
  void add(unsigned int _counter);
  void print(unsigned int _counter,const char * msg);
  void printAll(const char * msg);
private:
  myTimer(const myTimer &t); // this far not implemented - do not copy by accident
  unsigned int numCounters;
  unsigned long long int * counters;
  struct timespec last;
};


/***********************************
 ***   Template Implementation   ***
 **********************************/
template<typename T>
T mymin(T a, T b){
  return a>b?b:a;
}

template<typename T>
void printMatrixHost( T * host, int n, int m, const char * msg){
  std::cout << msg;
  int sym = mymin(m,n);
  int maxm=mymin(1000,sym);
  int maxn=mymin(1000,sym);
  if ( m<=0 || n<=0){
    return;
  }
  std::cout << "\nA=[";
  for (int j=0; j < maxm;j++){
    std::cout <<"[";
    std::cout << std::setw(12)<< host[n*j];
    for (int i =1; i<maxn;i++){
      std::cout << " , " <<std::setw(20) << std::setprecision(15)<< host[n*j+i];
      //printf(",\t%20.15e",host[j]);
    }
    //std::cerr << j <<" (" <<n<<","<<m<<") " << std::endl;
    if (j==n-1){
      std::cout << "]";
    } else {
      std::cout << "]\n   ";
    }
  }
  std::cout << "];\neig(A)\n";
}
template<typename T>
void printMatrixDev( T * data, int n, int m){
  T host[m*n];
  if ( m<=0 || n<=0){
    return;
  }
  else{
    cudaMemcpy( host, data, m*n*sizeof(*data), cudaMemcpyDeviceToHost );
    printMatrixHost(host,n,m,"");
  }
}


template<typename T>
void printMatrixDev( T * data, int m, int n, const char * msg){
  std::cout << msg;
  printMatrixDev(data,m,n);
}


template<typename T> void printDistances( T * r_h, T L, int N, T min, T max){
#ifndef DIM
  const int DIM=3;
#endif
  for (int i = 0;i<N;i++){
    for (int j=0; j < i;j++){
      T dr2=0;
      for (int k=0;k<DIM; k++){
	T tmp=r_h[DIM*i+k]-r_h[DIM*j+k];
	tmp-=L*rint(tmp/L);
	dr2+=tmp*tmp;
      }
      T drn=sqrt(dr2);
      if (min < drn && drn < max){
	std::cout << drn << "\n";
      }
    }
  }
}
template<typename T> void printDistances( T * r_h, T L, int N){
  printDistances(r_h, L, N,0,L);
}



#endif
