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



void ValueAdded(List All)
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



    Model.DataTransform();
    cout << "done " << "data update" << endl;
    Model.SigmaEst();
    cout << "done " << "sig" << endl;
    //-----------------------------

    Model.GetQH();
    cout << "done " << "qh" << endl;
    Model.Lambda();
    cout << "done " << "lmd" << endl;
    
    clock_t begin = clock();
    
    Model.Omega();
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    cout << elapsed_secs<<endl;
    
    cout << "done " << "omg" << endl;
    // Model.VA();
    // cout << "done " << "va" << endl;
    // mat Gamma = Model.Gamma;
    // Gamma = Gamma.t();
    // return Gamma;

  }








// [[Rcpp::export]]

void test1(){
  

  
  
  // std::vector<mat> l(25);
  // 
  // mat M(1000,1000, fill::randu);
  // 
  // 
  // int range = pow(10,2) ;
  // 
  // for(int i=0;i<range; i++){
  // 
  //   // cout<<i<<endl;
  // 
  //   l.push_back(M);
  // 
  // }
  // 
  // for(int i=0;i<range; i++){
  // 
  //   // cout<<i<<endl;
  // 
  //   mat  m = l[i];
  // 
  // }

// 
//   for(mat n:l){
// 
//     // cout<<i<<endl;
// 
//     m = n;
// 
//   }
  // int d = pow(10,3);
  // sp_mat m = speye<sp_mat>(d,d);
  // 
  // mat m2 = conv_to<mat>::from(m);
  // mat i = inv(m2);
  // 


  mat X = randu<mat>(2000,2000);
  
  clock_t begin = clock();
  
  
  mat inv_x = inv(X);
  

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  cout << elapsed_secs<<endl;
  //-----------------------------------------
    
  // clock_t begin2 = clock();
  // 
  // 
  // for(int i=0;i<range; i++){
  // 
  //   // cout<<i<<endl;
  //   mat EyeMat = eye<mat>(5*100,5*100);
  // 
  // }
  // 
  // clock_t end2 = clock();
  // double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
  // 
  // cout << elapsed_secs2<<endl;
  
}


// [[Rcpp::export]]

void test2(){

  List l;

  mat M(1000,1000, fill::randu);


  int range = pow(10,2) ;

  for(int i=0;i<range; i++){

   // cout<<i<<endl;

    l.push_back(M);

  }
  // 
  // 
  for(int i=0;i<range; i++){

    // cout<<i<<endl;

   mat  m = l[i];

  }
  
  // int d = pow(10,3);
  // mat m(d,d, fill::eye);
  // mat i = inv(m);
  // 

}
























