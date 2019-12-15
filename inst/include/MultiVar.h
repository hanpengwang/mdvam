#ifndef __MultiVarHeader__
#define __MultiVarHeader__
//#define ARMA_64BIT_WORD
//#include <Rcpp.h>
//#include <Rcpp.h>
#include "RcppArmadillo.h"
#include <algorithm>




// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

class MultiVar

  {

private:

//---------------Data&Params from original-------------
  int M; //number of Ys
  int N; //total obervations
  colvec Nj; //obervations for each category
  int K; //number of Xs (exclude constant regressor)
  int J; //number of Categories
  int DF; // degree of freedom
  List DataX, DataY;

//---------------Processed Data&Params -----------------
  //List ListZ;
  std::vector<mat> ListZ;
  mat Z;
  mat X;
  mat Y;
  mat WithinX;
  mat WithinY;
  mat MDiag;
  //Other Estimation;
  colvec Beta;
  colvec BetaGls;
  colvec e;
  double SigmaSquare;
  mat R;
  mat QH;
  mat LambdaTilde;
  mat LambdaUnarranged;
  sp_mat OmegaMat;
  //List OmegaList;
  std::vector<mat> OmegaList;

public:

  //-------------------------

  //MultiVar(List a, List b, int InputM, int InputN, colvec InputNj,
  //         int InputK, int InputJ, int InputDF);
  MultiVar();

  ~MultiVar();
  void SetData(List a, List b, int InputM, int InputN, colvec InputNj,
               int InputK, int InputJ, int InputDF);
  void UpdateData();
  void DataTransform();
  void SigmaEst();
  void GetQH();
  void Lambda();
  void Omega();
  void VA();
  mat Gamma;

  //some useful matrices

  mat Withinj(int& nj);
  mat Hj(int& nj);
  sp_mat Pjmj(int& m, int& nj);
  sp_mat Pjm(int& m, int& j, int& nj);
  mat bdiag(const mat& dmat, int &size);
  //mat bdiag(const List& ListMat);
  
  //sp_mat bdiag(const List& ListMat);
  sp_mat bdiag(const std::vector<mat>& ListMat);
  
  //some useful methods
  colvec SetDiff(colvec& x, colvec& y);
  mat FillTri(mat& FillMat, colvec& ValueVec, bool diag);

};



#endif
