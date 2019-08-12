#include "MultiVar.h"
#include <Rcpp.h>
#include "RcppArmadillo.h"

using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
  
MultiVar::MultiVar(List a, List b, List c)
  {
    DataY = a;
    DataX = b;
    DataJ = c; 
  }

MultiVar::~MultiVar(){}


void MultiVar::UpdateParameters(int InputM, int InputN, 
                                int InputK, int InputJ, int InputDF)
  
  { 
    M = InputM;
    N = InputN;
    K = InputK;
    J = InputJ;
    DF = InputDF;
  
  }

void MultiVar::DataTransform()
  {
    int FirstR = 0;
    int LastR;
    const int FirstC = 0;
    const int LastC_Z = FirstC + (K+1) * M - 1;
    const int LastC_X = FirstC + K * M - 1;
    const int LastC_Y = 0;
    
    for (int i=0; i<J; i++ )
      {
        int nj = Nj[i];
        mat wj = Withinj(nj);
        mat y = DataY[i]; y.set_size((nj * M), 1);
        mat withiny = wj * y;
        mat x = DataX[i]; x = MDiag * x;
        mat onevec(nj, 1, fill::ones);
        mat z = join_horiz(onevec, x); z = MDiag * z;
        mat withinx = wj * x;
        
        ListZ[i] = z;
        
        //fill each category matrices into variable matrices;
        
        LastR = FirstR + nj * M - 1;
        
        Z.submat(FirstR, FirstC, LastR, LastC_Z) = z;
        X.submat(FirstR, FirstC, LastR, LastC_X) = x;
        Y.submat(FirstR, FirstC, LastR, LastC_Y) = z;
        WithinX.submat(FirstR, FirstC, LastR, LastC_X) = withinx;
        WithinY.submat(FirstR, FirstC, LastR, LastC_Y) = withinx; 
        
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
            int nj2 = Nj[i1]; 
            mat Z2 = ListZ[i2];
            mat blocki = 0 - Z1 * (R * Z2.t());
            
            if(i2 == i1) blocki = eye<mat>(M * nj1, M*nj1);
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
        int FirstR_Gamma = j * M*(M+1)/2;
        int LastR_Gamma = (j + 1) * M*(M+1)/2 - 1;
        mat Alphaj = zeros<mat>(M*(M+1)/2,J);
        int AlphajCount = 0;
        int nj = Nj[j];
        int j2 = M * sum(Nj.subvec(0, j));
        if(j == 0) 
          {
          //int j1 =  M * as_scalar(sum(Nj.subvec(0, 0)));
          ej = e.subvec(0, j2);
          }
        else
          { 
          int j1 = M * as_scalar(sum(Nj.subvec(0, (j-1))));
          ej = e.subvec(j1, j2);
          }        

        mat Zj = ListZ[j];
      
        for (int p=1; p<=M; p++ )
          {
            sp_mat Pp = Pjm(p, j, nj);
            sp_mat Ppj = Pjmj(p, nj);
            colvec ejp = Ppj * ej;
            for (int q=p; q<=M; q++ )
              {
                sp_mat Pq = Pjm(q, j, nj);
                sp_mat Pqj = Pjmj(q, nj);
                mat QjTilde = eye<mat>(M*nj, M*nj) - Zj * (R * Zj.t());
                double vjpq = as_scalar(ejp.t() * Pqj * ej);
                double bjpq = SigmaSquare * sum((Pqj * (QjTilde * Pqj.t())).diag());
                mat Apq = (Pp * QH).t() * (Pq * QH);
                
                rowvec Alphajlpq(J);
                for (int d=0; d<=(J-1); d++ )
                  {
                   int d1 = d*M;
                   int d2 = d*M + M - 1;
                   mat Apqj = Apq.submat(d1,d1,d2,d2);
                   Alphajlpq(d) = Apqj(p-1,q-1);
                   
                  }
                v(vbCount) = vjpq;
                b(vbCount) = bjpq;
                vbCount++;
                Alphaj.row(AlphajCount) = Alphajlpq;
                AlphajCount++;
                
                
              }
          }
        for (int l=0; l<J; l++ )
          {
            
            int FirstC_Gamma = l * M*(M+1)/2;
            int LastC_Gamma = (l+1) * M*(M+1)/2 - 1;
            Gamma.submat(FirstR_Gamma, FirstC_Gamma, LastR_Gamma, LastC_Gamma) = diagmat(Alphaj.col(l));
          }
      }
    mat LambdaTilde_ = solve(Gamma,(v-b));
    
    LambdaTilde_.reshape((M*(M-1)/2 + M),J);
    
  
    //-------------------------------------rearrange lambdatilde---------------------------------------------------------
    
    int offset = 0;
    std::vector<int> Cols(M);
    Cols[0] = offset;
    int temp = 0;
    for(int i=1;i<=(M-1);i++)
      {
        
        offset = offset + M -temp;
        temp++;
        Cols[i] = offset;
      }
    
    std::vector<int> Cols2((M*(M+1)/2));
    std::iota(Cols2.begin(), Cols2.end(), 0);
    
    std::vector<int> diff;
    
    std::set_difference(Cols2.begin(), Cols2.end(), Cols.begin(), Cols.end(), 
                        std::inserter(diff, diff.begin()));
    
    std::vector<int> LambdaCols(Cols);
    LambdaCols.insert(LambdaCols.end(), diff.begin(), diff.end());
    
    LambdaTilde.set_size(LambdaTilde_.n_rows, LambdaTilde_.n_cols);
    
    for(int i; i<(M*(M+1)/2); i++)
      {
        LambdaTilde.col(i) = LambdaTilde_.col(LambdaCols[i]);
      }
    
  }


void MultiVar::Omega()
  {
  
    for(int i=0; i<J; i++)
      {
        int nj = Nj[i];
        colvec LambdaTildeJ = LambdaTilde.col(i);
        mat LambdaJ = zeros<mat>(M, M);
        
        int FirstI = 0;
        int LastI = M-1;
        for(int i2=0; i2<M; i2++ )
          {
            LambdaJ.submat(i2,i2,M-1,i2) = LambdaTildeJ.subvec(FirstI, LastI);
            FirstI = LastI + 1; 
            LastI = FirstI + (M-i2-2);
          }
        mat DiagLambdaJ = diagmat(LambdaJ);
        LambdaJ = LambdaJ + LambdaJ.t() - DiagLambdaJ;
        
        mat Jnj = ones<mat>(nj, nj);
        mat KronMat = kron(LambdaJ, Jnj);
        mat EyeMat = eye<mat>(M*nj,M*nj);
        mat Omegaj = KronMat + EyeMat*SigmaSquare;
        OmegaList[i] = inv(Omegaj);
        
      }
    OmegaMat = bdiag(OmegaList);
    
  }

void MultiVar::VA()
  {
    BetaGls = inv(Z.t() * OmegaMat * Z) * (Z.t() * OmegaMat * Y);
    colvec VGls = Y - Z * BetaGls;
    colvec Njm = Nj * M;
    int VGlsFirst = 0;
    int VGlsLast;
    for(int i=0;i<J;i++)
      {
        colvec TempLambdaJ = LambdaTilde.col(i);
        
        int VecLength = TempLambdaJ.n_rows; //15
        
          
        mat VCMLambdaJ = zeros<mat>(M,M); 
        int FirstI = 0;
        int LastI = FirstI + M - 1; 
        for(int i2=0; i2<M; i2++)
          {
            VCMLambdaJ.diag(i2) = TempLambdaJ.subvec(FirstI, LastI); 
            VCMLambdaJ.diag(-i2) = TempLambdaJ.subvec(FirstI, LastI);
            FirstI = LastI + 1; 
            LastI = FirstI + (M-i2-2);
            
             
            
          }
        int njm = Njm[i]; 
        int nj = Nj[i];
        VGlsLast = VGlsFirst + njm; 
        
        mat OmegaJ = OmegaList[i];
        mat Iotanj = ones<mat>(nj,nj);
        mat a = kron(VCMLambdaJ.t(), Iotanj);
        mat Gammaj = a * OmegaJ * VGls.subvec(VGlsFirst, VGlsLast);
        Gamma.col(i) = Gammaj;
        VGlsFirst = VGlsLast + 1; 
        
         
        
        
        
        
      }
  
    
    
  
  
  }










//------------------------some useful matrices------------------------------------------------------------------------------------------------------------------------------------

mat MultiVar::Withinj(int nj)
  {
    int MatSize = nj * M;
    mat Jnj = (ones<mat>(nj,nj)) / nj;
    mat wj = eye(MatSize, MatSize) - bdiag(Jnj, M) ;
    return wj;
  }


mat MultiVar::Hj(int nj)
  {
    mat OneMat = ones<mat>(1, nj);
    return bdiag(OneMat, M);
    
    
  }


sp_mat MultiVar::Pjmj(int &m, int &nj)
  {
    int n_row = nj;
    int n_col = (m - 1) * nj + nj + (M - m) * nj;
    sp_mat insert(n_row, n_col);
    insert.cols(((m - 1) * nj), ((m - 1) * nj + nj - 1)) = eye(nj, nj);
    
    return insert;
  }


sp_mat MultiVar::Pjm(int &m, int &j, int &nj)
  {
    int n_row = nj;
    int n_col = (N - nj) * M + (m - 1) * nj + nj + (M - m) * nj; 
    sp_mat p(n_row, n_col);
    sp_mat insert = Pjmj(m, nj);
    int FirstC = sum(Nj.subvec(0, (j - 1)));
    int LastC = FirstC + (m - 1) * nj + nj + (M - m) * nj - 1;
    p.cols(FirstC, LastC) = insert;
    
    return p;
  }


mat MultiVar::bdiag(const mat& dmat, int size)
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



mat MultiVar::bdiag(const List& ListMat)
{
  int size = ListMat.size();
  mat dmat = ListMat[0];
  
  mat bdm = zeros<mat>(dmat.n_rows * size, dmat.n_cols * size);
  for (int i = 0; i < size; ++i)
  {
    int FirstR = i * dmat.n_rows; 
    int FirstC = i * dmat.n_cols;
    int LastR = FirstR + dmat.n_rows - 1  ;
    int LastC = FirstC + dmat.n_cols - 1 ;
    mat BlockMat = ListMat[i];
    bdm.submat(FirstR, FirstC, LastR, LastC) = BlockMat;
  }
  
  return bdm;
} 


