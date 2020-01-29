//#include <Rcpp.h>
//#define ARMA_64BIT_WORD

#include "RcppArmadillo.h"
#include "MultiVar.h"
#include <math.h> 
#include <ctime>

using namespace Rcpp;
using namespace arma;
using namespace std;



// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]



arma::mat ValueAdded(List All)
  {

    List DataY = All[0];
    List DataX = All[1];
    int M = All[2];
    int N = All[3];
    arma::colvec Nj = All[4];
    int K = All[5];
    int J = All[6];
    int DF = All[7];

    //sp_mat SpOmega= All[8];

    MultiVar Model;

    Model.SetData(DataY, DataX, M, N, Nj,
                  K, J, DF);
    Model.UpdateData();

    clock_t begin = clock();

    Model.DataTransform();
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    //cout << elapsed_secs<<endl;
    
    
    
    //cout << "done " << "data update" << endl;
    Model.SigmaEst();
    //cout << "done " << "sig" << endl;
    //-----------------------------

    Model.GetQH();
    //cout << "done " << "qh" << endl;
    Model.Lambda();
    //cout << "done " << "lmd" << endl;
    
    // clock_t begin = clock();
    
    Model.Omega();
    
    // clock_t end = clock();
    // double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    // cout << elapsed_secs<<endl;

    
    //cout << "done " << "omg" << endl;
    Model.VA();
    //cout << "done " << "va" << endl;
    mat Gamma = Model.Gamma;
    Gamma = Gamma.t();
    
    return Gamma;
    
  }










// void test3(){
//   mat A(1000, 50, fill::randu);
//   mat C;
//   
//   clock_t begin = clock();
//   
//   for(int i = 0; i < 1000; i++){
//     C = A * A.t();
//   }
//   
//   clock_t end = clock();
//   double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//   
//   cout << elapsed_secs<<endl;
//   
//   
//   
//   clock_t begin2 = clock();
// 
//   std::vector<mat> ListC;
//   ListC.reserve(20);
// 
//   for(int i = 0; i < 100; i++){
//       ListC.push_back(A * A.t());
//   }
// 
//   for(int i = 0; i < 1000; i++){
//       C = ListC[(21+i)%20];
//   }
//   
//   clock_t end2 = clock();
//   double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
//   
//   cout << elapsed_secs2<<endl;
// 
// 
// }




















