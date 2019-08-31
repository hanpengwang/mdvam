//#include <Rcpp.h>
#include "RcppArmadillo.h"
#include "MultiVar.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
List testfunc(List x) {
  //x.push_back(1);
  List y;
  y.push_back(x[0]);
  return y;


}




/*

// [[Rcpp::export]]
void printsomething(){

  Rcout << "hello world/n";
  testmodel test;
  test.Printyou();

}

*/



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
    return Model.Gamma;



  }
