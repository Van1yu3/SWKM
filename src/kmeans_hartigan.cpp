
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(cpp11)]]
// Calculate euclidean distance (square)
mat CalEucDist(mat X,mat C){
  mat Dis(X.n_rows,C.n_rows); //n by p
  mat temp(size(X));
  for(unsigned int i=0;i<C.n_rows;++i){
    temp=X.each_row()-C.row(i);
    temp=temp%temp;
    Dis.col(i)=sum(temp,1);
  }
  return Dis;
}
// Calculate nearest and 2nd nearest cluster
void CalCl(mat Dis,ivec& cl,ivec& cl2){
  rowvec dis_rowi;
  for(unsigned int i=0;i<Dis.n_rows;++i){
    dis_rowi=Dis.row(i);
    if(dis_rowi[1]<dis_rowi[0]){
      cl[i]=1;
      cl2[i]=0;
    }else{
      cl[i]=0;
      cl2[i]=1;
    }
    for(unsigned int j=2;j<Dis.n_cols;++j){
      if(dis_rowi[j]<dis_rowi[cl[i]]){
        cl2[i]=cl[i];
        cl[i]=j;
      }else if(dis_rowi[j]<dis_rowi[cl2[i]]){
        cl2[i]=j;
      }
    }
  }
  ivec check=unique(cl);
  if(check.n_elem!=Dis.n_cols)
    Rf_error("empty cluster: try a better set of initial centers.");
}

// [[Rcpp::export]]
List kmeans_hartigan(NumericMatrix data, NumericMatrix centers, NumericVector weight){
  mat X=as<mat>(data);
  mat C=as<mat>(centers);
  vec W=as<vec>(weight);
  int K=centers.nrow();
  int n=data.nrow();
  int niter=0;
  int maxiter=10;
  ivec cl(n);
  ivec cl2(n);
  ivec live;
  live.ones(K);
  double wcss=0;

  // Step 1-3: initialization, update cluster centers and record the number of points within each cluster
  mat Dis=CalEucDist(X,C);
  CalCl(Dis,cl,cl2);
  vec Wsum=zeros(K);
  ivec num_cl;
  num_cl.zeros(K);
  int index,index2;
  C.zeros();
  for(int i=0;i<n;++i){
    index=cl[i];
    C.row(index) += W[i] * X.row(i);
    Wsum[index] += W[i];
    ++(num_cl[index]);
  }
  for(int i=0;i<K;++i){
    C.row(i) /= Wsum[i];
  }

  // Step 4: optimal-transfer stage
  double R1;
  double R2;
  double R2_tmp;
  int R2_index=0;
  vec ct_dif;
  while(niter<maxiter){
    Dis=CalEucDist(X,C);
    mat C_old(C);
    for(int i=0;i<n;++i){
      index=cl[i];
      R2=DBL_MAX;
      if(num_cl[index]>1){
        R1=Dis(i,index)*Wsum[index]/(Wsum[index]-W[i]);
        if(!std::isfinite(R1))
          Rf_error("some weights are so small that they are ignored. Please try better candidate weights.");
        if(live[index]==1){
          for(int nc=0;nc<K;++nc){
            if(nc!=index){
              R2_tmp=Dis(i,nc)*Wsum[nc]/(Wsum[nc]+W[i]);
              if(R2_tmp<R2) {
                R2=R2_tmp;
                R2_index=nc;
              }
            }
          }
        }else{
          for(int nc=0;nc<K;++nc){
            if(live[nc]==1 && nc!=index){
              R2_tmp=Dis(i,nc)*Wsum[nc]/(Wsum[nc]+W[i]);
              if(R2_tmp<R2) {
                R2=R2_tmp;
                R2_index=nc;
              }
            }
          }
        }
        if(R1<=R2) cl2[i]=R2_index;
        else{
          cl2[i]=cl[i];
          cl[i]=R2_index;

          C.row(R2_index) *= Wsum[R2_index];
          C.row(cl2[i]) *= Wsum[cl2[i]];

          Wsum[R2_index] += W[i];
          Wsum[cl2[i]] -= W[i];

          C.row(R2_index) += W[i]*X.row(i);
          C.row(cl2[i]) -= W[i] * X.row(i);

          C.row(R2_index) /= Wsum[R2_index];
          C.row(cl2[i]) /= Wsum[cl2[i]];

          live[R2_index]=1;
          live[cl2[i]]=1;

          ++(num_cl[R2_index]);
          --(num_cl[cl2[i]]);

          Dis=CalEucDist(X,C);
        }
      }
    }
    // Step 5: stop if the live set is empty
    ct_dif=sum(abs(C-C_old),1);
    for(int nc=0;nc<K;++nc){
      if(ct_dif[nc]==0) live[nc]=0;
    }
    if(!any(live)) break;
    // Step 6: quick-transfer stage
    int QTRAN=1;
    while(QTRAN){
      QTRAN=0;
      for(int i=0;i<n;++i){
        index=cl[i];
        index2=cl2[i];
        if((live[index]==1 || live[index2]==1) && num_cl[index]>1){
          R1=Dis(i,index)*Wsum[index]/(Wsum[index]-W[i]);
          R2=Dis(i,index2)*Wsum[index2]/(Wsum[index2]+W[i]);
          if(R1>=R2){
            cl[i]=index2;
            cl2[i]=index;

            C.row(index2) *= Wsum[index2];
            C.row(index) *= Wsum[index];

            Wsum[index2] += W[i];
            Wsum[index] -= W[i];

            C.row(index2) += W[i]*X.row(i);
            C.row(index) -= W[i] * X.row(i);

            C.row(index2) /= Wsum[index2];
            C.row(index) /= Wsum[index];

            live[index2]=1;
            live[index]=1;

            ++(num_cl[index2]);
            --(num_cl[index]);

            Dis=CalEucDist(X,C);
            QTRAN=1;
          }
        }
      }
    }
    ++niter;
  }
  // Compute WCSS
  for(int i=0;i<n;++i){
    index=cl[i];
    wcss += Dis(i,index) * W[i] ;
  }
  wcss /= sum(weight);
  ++cl;
  return List::create(Named("centers")=wrap(C),
                      Named("cluster")=wrap(cl),
                      Named("weight")=weight,
                      Named("wcss")=wrap(wcss));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
set.seed(1)
require(mvtnorm)
  n <- 60  #sample size
p <- 1000 #dimension of features
q <- 50  #dimension of cluster-specific features
mu <- 0.8
MU <- c(0,-mu,mu)
data <- rbind(rmvnorm(n/3,rep(0,p)),rmvnorm(n/3,c(rep(-mu,q),rep(0,p-q))),
              rmvnorm(n/3,c(rep(mu,q),rep(0,p-q))))
  set.seed(1)
  centers <- data[sample(n,3),]
  res2 <- kmeans_hartigan(data,centers,weight=(rep(1,n)))
  set.seed(1)
  res <- kmeans.weight(data,3,rep(1,n),centers)
*/
