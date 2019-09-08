//#include <Rcpp.h>
//#define ARMA_64BIT_WORD

#include "RcppArmadillo.h"
#include "MultiVar.h"



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



    Model.DataTransform();
    cout << "done " << "data update" << endl;
    Model.SigmaEst();
    cout << "done " << "sig" << endl;
    //-----------------------------

    Model.GetQH();
    cout << "done " << "qh" << endl;
    Model.Lambda();
    cout << "done " << "lmd" << endl;
    Model.Omega();
    cout << "done " << "omg" << endl;
    Model.VA();
    cout << "done " << "va" << endl;
    mat Gamma = Model.Gamma;
    Gamma = Gamma.t();
    return Gamma;

  }






