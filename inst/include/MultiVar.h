#ifndef __MultiVarHeader__
#define __MultivarHeader__

//#include <Rcpp.h>
#include <Rcpp.h>
#include "RcppArmadillo.h"


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


class MultiVar
  
  {
  
private:
  
  int M; //number of Ys 
  int N; //total obervations
  colvec Nj; //obervations for each category
  int K; //number of Xs (exclude constant regressor)
  int J; //number of Categories
  int DF; // degree of freedom 
  colvec Beta; 
  colvec BetaGls;
  colvec e;
  double SigmaSquare;
  List ListZ = no_init(J) ;
  mat Z = zeros<mat>((N * M), ((K+1)*M));
  mat X = zeros<mat>((N * M), (K*M));
  mat Y = zeros<mat>((N * M), 1);
  mat WithinX = zeros<mat>((N * M), (K*M));
  mat WithinY = zeros<mat>((N * M), 1);
  List DataX, DataY, DataJ; 
  mat MDiag = eye<mat>(M,M);
  //Other Estimation Matrices;
  mat R;
  mat QH = zeros<mat>((N * M), Nj.n_rows * M);
  mat LambdaTilde;
  mat OmegaMat;
  List OmegaList;
  mat Gamma = zeros<mat>(J, M);
public:
  
  MultiVar(List a, List b, List c);
  ~MultiVar();
  void UpdateParameters(int InputM, int InputN, 
                        int InputK, int InputJ, int InputDF);
  void DataTransform();
  void SigmaEst();
  void GetQH();
  void Lambda();
  void Omega();
  void VA();
  
  
  //some useful matrices
  
  mat Withinj(int nj);
  mat Hj(int nj); 
  sp_mat Pjmj(int &m, int &nj);
  sp_mat Pjm(int &m, int &j, int &nj);
  mat bdiag(const mat& dmat, int size);
  mat bdiag(const List& ListMat);



};






























#endif
