//#define ARMA_64BIT_WORD
#include "MultiVar.h"
#include "RcppArmadillo.h"
#include <algorithm>

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]


MultiVar::MultiVar(){}

MultiVar::~MultiVar(){}

void MultiVar::SetData(List a, List b, int InputM, int InputN, colvec InputNj,
                       int InputK, int InputJ, int InputDF)
  {
    DataY = a;
    DataX = b;
    M = InputM; //number of Ys
    N = InputN; //total obervations
    Nj = InputNj; //obervations for each category
    K = InputK; //number of Xs (exclude constant regressor)
    J = InputJ; //number of Categories
    DF = InputDF; // degree of freedom

  }


void MultiVar::UpdateData()

  {
    Z = zeros<mat>((N * M), ((K+1)*M));
    X = zeros<mat>((N * M), (K*M));
    Y = zeros<mat>((N * M), 1);
    WithinX = zeros<mat>((N * M), (K*M));
    WithinY = zeros<mat>((N * M), 1);
    MDiag = eye<mat>(M,M); //* later check for opt
    Gamma = zeros<mat>(M, J); //* later check for opt
    QH = zeros<mat>((N * M), Nj.n_rows * M);  //* later check for opt

  }

void MultiVar::DataTransform()
  {
    int FirstR = 0;
    int LastR;
    const int FirstC = 0;
    const int LastC_Z = FirstC + (K+1) * M - 1;
    const int LastC_X = FirstC + K * M - 1;
    const int LastC_Y = 0;
    
    ListZ.reserve(J); //## initialize size of ListZ; 
    
    for (int i=0; i<J; i++ )
      {
        int nj = Nj[i];
        mat wj = Withinj(nj);
        mat y = DataY[i];
        y.set_size((nj * M), 1);
        mat withiny = wj * y;
        mat x = DataX[i];
        mat onevec(nj, 1, fill::ones);
        mat z = join_horiz(onevec, x);
        //#
        //## x = kron(MDiag, x);
        //## z = kron(MDiag, z);
        x = bdiag(x, M);
        z = bdiag(z, M);
        //#
        mat withinx = wj * x;
        ListZ.push_back(z);



        //fill each category matrices into variable matrices;

        LastR = FirstR + nj * M - 1;

        Z.submat(FirstR, FirstC, LastR, LastC_Z) = z;
        //cout << "done" << "z" << endl;
        X.submat(FirstR, FirstC, LastR, LastC_X) = x;
        //cout << "done" << "x" << endl;
        Y.submat(FirstR, FirstC, LastR, LastC_Y) = y;
        //cout << "done" << "y" << endl;
        WithinX.submat(FirstR, FirstC, LastR, LastC_X) = withinx;
        //cout << "done" << "wx" << endl;
        WithinY.submat(FirstR, FirstC, LastR, LastC_Y) = withiny;
        //cout << "done" << "wy" << endl;

        FirstR = LastR + 1;


    }



  }


void MultiVar::SigmaEst()
  {

    colvec WithinBeta = solve((WithinX.t() * WithinX) , (WithinX.t() * WithinY));
    colvec res = Y - X * WithinBeta;
    colvec wres = WithinY - WithinX * WithinBeta;
    mat CrossProd = res.t() * wres;
    double rss = as_scalar(CrossProd);
    SigmaSquare = rss / DF;
    Beta = solve((Z.t() * Z) , (Z.t() * Y));
    e = Y - Z * Beta;
    R = inv((Z.t() * Z));

}

void MultiVar::GetQH()
  {
    int FirstR = 0;
    int FirstC;
    int LastR;
    int LastC;
    for (int i1=0; i1<J; i1++ )
      {

        int nj1 = Nj[i1];
        mat Z1 = ListZ[i1];
        LastR = FirstR + nj1 * M - 1;

        for (int i2=0; i2<J; i2++ )
          {
            int nj2 = Nj[i2];
            mat Z2 = ListZ[i2];
            mat blocki = 0 - Z1 * (R * Z2.t()); 
            if(i2 == i1) blocki = eye<mat>(M * nj1, M*nj1) + blocki;
            mat QHj1j2 = blocki * Hj(nj2).t();
            FirstC = i2 * M;
            LastC = FirstC + M -1;
            QH.submat(FirstR, FirstC, LastR, LastC) = QHj1j2;

        }
        FirstR = LastR + 1;
    }

}

void MultiVar::Lambda()
  {
    int vbCount = 0;
    colvec v = zeros<colvec>((J*M*(M+1)/2), 1);
    colvec b = zeros<colvec>((J*M*(M+1)/2), 1);
    mat Gamma = zeros<mat>((J*M*(M+1)/2), (J*M*(M+1)/2));
    colvec ej;

    for (int j=0; j<J; j++ )
      {
        int AlphajCount = 0;
        int FirstR_Gamma = j * M*(M+1)/2;
        int LastR_Gamma = (j + 1) * M*(M+1)/2 - 1;
        mat Alphaj = zeros<mat>(M*(M+1)/2,J);
        int nj = Nj[j];
        int j2 = M * sum(Nj.subvec(0, j));

        //cout<<"j2 is " << j2 <<endl;


        if(j == 0)
          {
          //int j1 =  M * as_scalar(sum(Nj.subvec(0, 0)));
          ej = e.subvec(0, (j2-1));

          //cout<<"ej size " << ej.size() <<endl;

          }

        else
          {

          int j1 = M * sum(Nj.subvec(0, j-1));

          //cout<<"j1 " << j1 <<endl;

          ej = e.subvec(j1, j2-1);

          //cout<<"ej size " << ej.size() <<endl;

          }

        mat Zj = ListZ[j];

        for (int p=1; p<=M; p++ )
          {
            //cout<<"going here 222--------- " <<endl;

            sp_mat Pp = Pjm(p, j, nj);
            sp_mat Ppj = Pjmj(p, nj);
            colvec ejp = Ppj * ej;

            //cout<<"finish p " << p <<endl;


            for (int q=p; q<=M; q++ )
              {
                sp_mat Pq = Pjm(q, j, nj);
                sp_mat Pqj = Pjmj(q, nj);
                mat QjTilde = eye<mat>(M*nj, M*nj) - Zj * (R * Zj.t());  //* later check for opt
                double vjpq = as_scalar(ejp.t() * Pqj * ej);
                double bjpq = SigmaSquare * sum((Pqj * (QjTilde * Ppj.t())).diag());
                mat Apq = (Pp * QH).t() * (Pq * QH);

                //Alpha.push_back(Apq);

                rowvec Alphajlpq(J);

                //cout<<"finish q " << q <<endl;

                for (int d=0; d<=(J-1); d++ )
                  {
                   int d1 = d*M;
                   int d2 = d*M + M - 1;
                   mat Apqj = Apq.submat(d1,d1,d2,d2);
                   Alphajlpq(d) = Apqj(p-1,q-1);
                   //cout << "--------------------"<<endl;
                   //cout<< Apqj<< endl;

                   //cout<< d1 << " " << d2 << endl;

                   //cout<<"finish d " << d <<endl;

                  }
                v(vbCount) = vjpq;
                b(vbCount) = bjpq;
                vbCount++;
                Alphaj.row(AlphajCount) = Alphajlpq;
                AlphajCount++;


              }
            //Alpha.push_back(Alphaj);
          }
        for (int l=0; l<J; l++ )
          {

            int FirstC_Gamma = l * M*(M+1)/2;
            int LastC_Gamma = (l+1) * M*(M+1)/2 - 1;
            Gamma.submat(FirstR_Gamma, FirstC_Gamma, LastR_Gamma, LastC_Gamma) = diagmat(Alphaj.col(l));
            //mat DiagAlpha(Alphaj.col(l).n_rows,Alphaj.col(l).n_rows);  DiagAlpha.eye();
            //Gamma.submat(FirstR_Gamma, FirstC_Gamma, LastR_Gamma, LastC_Gamma) = Alphaj.col(l)[0] * DiagAlpha; Alpha.push_back(Alphaj.col(l)[0] * DiagAlpha);
          }

        //cout<<"done " << j <<endl;

      }


    mat LambdaTilde_ = solve(Gamma, (v-b));


    LambdaTilde_.reshape((M*(M-1)/2 + M),J);

    //-------------------------------------rearrange lambdatilde---------------------------------------------------------

    int offset = 0;
    colvec Rows(M);
    Rows(0) = offset;
    int temp = 0;
    for(int i=1;i<=(M-1);i++)
      {

        offset = offset + M -temp;
        temp++;
        Rows(i) = offset;


      }

    colvec sequence = linspace(0, LambdaTilde_.n_rows - 1, LambdaTilde_.n_rows);
    colvec difference = SetDiff(sequence, Rows);
    colvec AllRows = join_cols(Rows, difference);


    //LambdaTilde = LambdaTilde_.rows(AllCols);
    /*
    std::vector<int> Cols2((M*(M+1)/2));
    std::iota(Cols2.begin(), Cols2.end(), 0);

    std::vector<int> diff;

    std::set_difference(Cols2.begin(), Cols2.end(), Cols.begin(), Cols.end(),
                        std::inserter(diff, diff.begin()));

    std::vector<int> LambdaCols(Cols);
    LambdaCols.insert(LambdaCols.end(), diff.begin(), diff.end());

    LambdaTilde.set_size(LambdaTilde_.n_rows, LambdaTilde_.n_cols);
*/
    LambdaTilde.set_size((M*(M+1)/2), J);
    for(int i=0; i<(M*(M+1)/2); i++)
      {
        LambdaTilde.row(i) = LambdaTilde_.row(AllRows[i]);
      }
    LambdaUnarranged = LambdaTilde_;

  }


void MultiVar::Omega()
  {
    mat LambdaJ_template = zeros<mat>(M, M);
    for(int i=0; i<J; i++)
      {
        //cout << "going here 1--------------------" << endl;
        int nj = Nj[i];
        colvec LambdaTildeJ = LambdaUnarranged.col(i);
        mat LambdaJ = LambdaJ_template;
        //mat LambdaJ = zeros<mat>(M, M);
        LambdaJ = FillTri(LambdaJ, LambdaTildeJ, true);

        //cout << "going here 2--------------------" << endl;

        /*
        int FirstI = 0;
        int LastI = M-1;
        cout << "going here 2--------------------" << endl;
        for(int i2=0; i2<M; i2++ )
          {
            LambdaJ.submat(i2,i2,M-1,i2) = LambdaTildeJ.subvec(FirstI, LastI);
            FirstI = LastI + 1;
            LastI = FirstI + (M-i2-2);
            //cout << "going here 3--------------------" << endl;
          }
         */
        mat DiagLambdaJ = diagmat(LambdaJ);
        LambdaJ = LambdaJ + LambdaJ.t() - DiagLambdaJ;
        //Lambdaj.push_back(LambdaJ);

        //cout << "going here 4--------------------" << endl;

        mat Jnj = ones<mat>(nj, nj);
        //cout << "going here 5--------------------" << endl;
        mat KronMat = kron(LambdaJ, Jnj);
        //cout << "going here 6--------------------" << endl;
        mat EyeMat = eye<mat>(M*nj,M*nj);
        mat Omegaj = KronMat + EyeMat*SigmaSquare;
        //cout << "going here 7--------------------" << endl;
        //clock_t begin = clock();
        
        OmegaList.push_back(inv(Omegaj));
        //clock_t end = clock();
        //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        
        //cout << elapsed_secs<<endl;
        //cout << "going here 8--------------------" << endl;

      }

    //mat omega0 = OmegaList[0];
    //cout << omega0.n_rows<< " " << omega0.n_cols<< endl;

    OmegaMat = bdiag(OmegaList);
    //Omegabdiag(OmegaList);

  }

void MultiVar::VA()
  {
    BetaGls = inv(Z.t() * OmegaMat * Z) * (Z.t() * OmegaMat * Y);
    colvec VGls = Y - Z * BetaGls;
    colvec Njm = Nj * M;
    int VGlsFirst = 0;
    int VGlsLast;
    //cout << "going here 3--------------------" << endl;
    for(int i=0;i<J;i++)
      {
        colvec TempLambdaJ = LambdaUnarranged.col(i);

        int VecLength = TempLambdaJ.n_elem; //15


        mat VCMLambdaJ = zeros<mat>(M,M);
        int FirstI = 0;
        int LastI = FirstI + M - 1;
        for(int i2=0; i2<M; i2++)
          {
            VCMLambdaJ.submat(i2,i2,M-1,i2) = TempLambdaJ.subvec(FirstI, LastI);
            VCMLambdaJ.submat(i2,i2,i2,M-1) = TempLambdaJ.subvec(FirstI, LastI).t();
            //VCMLambdaJ.diag(-i2) = TempLambdaJ.subvec(FirstI, LastI);
            FirstI = LastI + 1;
            LastI = FirstI + (M-i2-2);



          }
        //Lambdaj.push_back(VCMLambdaJ);

        int njm = Njm[i];
        int nj = Nj[i];
        VGlsLast = VGlsFirst + njm - 1;

        mat OmegaJ = OmegaList[i];
        mat Iotanj = ones<rowvec>(nj);
        //cout << VCMLambdaJ.n_rows<< " " << VCMLambdaJ.n_cols  << endl;
        mat a = kron(VCMLambdaJ.t(), Iotanj);
        //cout << a.n_rows<< " " << a.n_cols  << endl;
        mat Gammaj = a * OmegaJ * VGls.subvec(VGlsFirst, VGlsLast);
        //cout << "going here 2--------------------" << endl;
        Gamma.col(i) = Gammaj;
        VGlsFirst = VGlsLast + 1;

        //cout << "done " << i<< endl;




      }





  }










//------------------------some useful matrices------------------------------------------------------------------------------------------------------------------------------------

mat MultiVar::Withinj(int& nj)
  {
    //## int MatSize = nj * M;
    //## mat Jnj = (ones<mat>(nj,nj)) / nj;
    //##mat wj = eye(MatSize, MatSize) - bdiag(Jnj, M) ;
    mat Jnj(nj, nj); 
    double nj_1 = - 1 / nj;
    Jnj.fill(nj_1);
    Jnj.diag() += 1; 
    mat wj = bdiag(Jnj, M);
    //##----
    return wj;
  }

  
mat MultiVar::Hj(int& nj)
  {
    mat OneMat = ones<mat>(1, nj);
    return bdiag(OneMat, M);
  }


sp_mat MultiVar::Pjmj(int& m, int& nj)
  {
    int n_row = nj;
    int n_col = (m - 1) * nj + nj + (M - m) * nj;
    sp_mat insert(n_row, n_col);
    insert.cols(((m - 1) * nj), ((m - 1) * nj + nj - 1)) = eye(nj, nj);

    return insert;
  }


sp_mat MultiVar::Pjm(int& m, int& j, int& nj)
  {
    int FirstC;
    int n_row = nj;
    int n_col = (N - nj) * M + (m - 1) * nj + nj + (M - m) * nj;
    sp_mat p(n_row, n_col);
    sp_mat insert = Pjmj(m, nj);
    if (j == 0)
      {
        FirstC = 0;
      }
    else
      {
        FirstC = M * sum(Nj.subvec(0, (j - 1)));
      }
    int LastC = FirstC + (m - 1) * nj + nj + (M - m) * nj - 1;
    p.cols(FirstC, LastC) = insert;

    return p;
  }


mat MultiVar::bdiag(const mat& dmat, int& size)
  {
    mat bdm = zeros<mat>(dmat.n_rows * size, dmat.n_cols * size);
    for (int i = 0; i < size; ++i)
    {
      int FirstR = i * dmat.n_rows;
      int FirstC = i * dmat.n_cols;
      int LastR = FirstR + dmat.n_rows - 1  ;
      int LastC = FirstC + dmat.n_cols - 1 ;
      bdm.submat(FirstR, FirstC, LastR, LastC) = dmat;
    }

    return bdm;
  }



//sp_mat MultiVar::bdiag(const List& ListMat)
sp_mat MultiVar::bdiag(const std::vector<mat>& ListMat)
{
  int FirstRC = 0;
  int LastRC;
  colvec Mnj = Nj * M;
  sp_mat bdm = sp_mat(M*N, M*N);
  for (int i = 0; i < J; ++i)
  {
    LastRC = FirstRC + Mnj[i] - 1;
    mat Omegaj = ListMat[i];
    bdm.submat(FirstRC, FirstRC, LastRC, LastRC) = Omegaj;
    FirstRC = LastRC + 1;
  }

  return bdm;
}

arma::colvec MultiVar::SetDiff(arma::colvec& x, arma::colvec& y) {

  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;

  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));

  return arma::conv_to<arma::colvec>::from(out);
}


mat MultiVar::FillTri(mat& FillMat, colvec& ValueVec, bool diag)
  {
     /*
      * FillMat is a symmetric matrix;
      */
    //int nrow = FillMat.n_rows;
    int ncol = FillMat.n_cols;
    
    //int VecLen = ValueVec.n_elem;
    //if(nrow != ncol) stop("matrix must be symetric");
    //if(diag == true && (nrow**2 - 1/2 * nrow)!= VecLen  ) stop("lengths are different")
    //if(diag == false && (nrow**2)!= VecLen  ) stop("lengths are different")

    // vertical fill , fill columm 1 then 2 ,3 ...........
    int VecCount = 0;
    // if(diag == true )
    //   {
        for(uword pcol=0; pcol<ncol; pcol++)
          {
            for(uword prow=pcol; prow<ncol; prow++)
              {
                FillMat(prow, pcol) = ValueVec[VecCount];
                VecCount++;
              }
          }
        
        /*
         * alternative way to fillmat
        
        for(uword pcol=0; pcol<ncol; pcol++)
        {
          uword indexer = pcol * 5;
          for(uword prow=pcol; prow<ncol; prow++)
          {
            uword innerIndexer = indexer + prow ;
            FillMat.at(innerIndexer) = ValueVec[VecCount];
            VecCount++;
          }
        }
      
         */
        
        // }

    return FillMat;



  }

