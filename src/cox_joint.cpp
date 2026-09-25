#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_joint_lambda(double t) {
  double res = 0.5*t;
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_joint_Lambda(double t) {
  double res = 0.25*pow(t,2);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_joint_log_normal_pdf(double t, double m, double s, double pi) {
  double res = -pow((t-m)/s,2)/2-log(2*pi)/2-log(s);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_joint_d_log_normal_pdf(double t, double m, double s) {
  double res = -(t-m)/pow(s,2);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_joint_dd_log_normal_pdf(double s) {
  double res = -1/pow(s,2);
  return res;
}

double cpp_cox_joint_dm_log_normal_pdf(double t, double m, double s) {
  double res = (t-m)/pow(s,2);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_joint_ddm_log_normal_pdf(double s) {
  double res = -1/pow(s,2);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_joint_est_Lambda(double t, arma::mat data_compl, arma::colvec est_compl, arma::colvec ipw_compl) {
  int p = est_compl.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec Y = data_compl.col(1);
  arma::colvec D = data_compl.col(2);
  arma::colvec X = data_compl.col(3);
  arma::mat Z = data_compl.submat(0,5,m-1,p+3);
  
  double res = 0;
  for(int i=0; i<m; i++) {
    if((Y[i]<=t)&&(D[i]==1)) {
      double temp = 0;
      for(int j=0; j<m; j++) {
        if(Y[j]>=Y[i]) {
          temp += ipw_compl[j]*exp(est_compl[0]*X[j]+as_scalar(Z.row(j)*est_compl.subvec(1,p-1)));
        }
      }
      res += ipw_compl[i]/temp;
    }
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_joint_est_lambda(arma::mat data_compl, arma::colvec est_compl, arma::colvec ipw_compl) {
  int p = est_compl.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec Y = data_compl.col(1);
  arma::colvec D = data_compl.col(2);
  arma::colvec X = data_compl.col(3);
  arma::mat Z = data_compl.submat(0,5,m-1,p+3);
  
  arma::colvec res = zeros(m);
  for(int i=0; i<m; i++) {
    if(D[i]==1) {
      double temp = 0;
      for(int j=0; j<m; j++) {
        if(Y[j]>=Y[i]) {
          temp += ipw_compl[j]*exp(est_compl[0]*X[j]+as_scalar(Z.row(j)*est_compl.subvec(1,p-1)));
        }
      }
      res[i] = ipw_compl[i]/temp;
    }
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_joint_s0_true(arma::colvec beta, arma::colvec t, arma::mat data_compl, arma::colvec est_compl, arma::colvec ipw_compl) {
  int p = beta.n_elem;
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec X_compl = data_compl.col(3);
  arma::mat Z_compl = data_compl.submat(0,5,m-1,p+3);
  
  double sipw = sum(ipw_compl);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double Lambda_ti = cpp_cox_joint_est_Lambda(t[i],data_compl,est_compl,ipw_compl);
    double temp = 0;
    for(int j=0; j<m; j++) {
      double XZbj = beta[0]*X_compl[j]+as_scalar(Z_compl.row(j)*beta.subvec(1,p-1));
      temp += ipw_compl[j]*exp(XZbj-Lambda_ti*exp(XZbj));
    }
    res[i] = temp/sipw;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_s1_true(arma::colvec beta, arma::colvec t, arma::mat data_compl, arma::colvec est_compl, arma::colvec ipw_compl) {
  int p = beta.n_elem;
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec X_compl = data_compl.col(3);
  arma::mat Z_compl = data_compl.submat(0,5,m-1,p+3);
  
  double sipw = sum(ipw_compl);
  
  arma::mat res = zeros(n,p);
  for(int i=0; i<n; i++) {
    double Lambda_ti = cpp_cox_joint_est_Lambda(t[i],data_compl,est_compl,ipw_compl);
    double tempX = 0;
    arma::rowvec tempZ = zeros<arma::rowvec>(p-1);
    for(int j=0; j<m; j++) {
      double XZbj = beta[0]*X_compl[j]+as_scalar(Z_compl.row(j)*beta.subvec(1,p-1));
      tempX += ipw_compl[j]*X_compl[j]*exp(XZbj-Lambda_ti*exp(XZbj));
      tempZ += ipw_compl[j]*Z_compl.row(j)*exp(XZbj-Lambda_ti*exp(XZbj));
    }
    res(i,0) = tempX/sipw;
    res.submat(i,1,i,p-1) = tempZ/sipw;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_joint_s0_org(arma::colvec beta, arma::colvec t, arma::mat data_compl, arma::colvec ipw_compl) {
  int p = beta.n_elem;
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec Y_compl = data_compl.col(1);
  arma::colvec X_compl = data_compl.col(3);
  arma::mat Z_compl = data_compl.submat(0,5,m-1,p+3);
  
  double sipw = sum(ipw_compl);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp = 0;
    for(int j=0; j<m; j++) {
      if(Y_compl[j]>=t[i]){
        temp += ipw_compl[j]*exp(beta[0]*X_compl[j]+as_scalar(Z_compl.row(j)*beta.subvec(1,p-1)));
      }
    }
    res[i] = temp/sipw;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_s1_org(arma::colvec beta, arma::colvec t, arma::mat data_compl, arma::colvec ipw_compl) {
  int p = beta.n_elem;
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec Y_compl = data_compl.col(1);
  arma::colvec X_compl = data_compl.col(3);
  arma::mat Z_compl = data_compl.submat(0,5,m-1,p+3);
  
  double sipw = sum(ipw_compl);
  
  arma::mat res = zeros(n,p);
  for(int i=0; i<n; i++) {
    double tempX = 0;
    arma::rowvec tempZ = zeros<arma::rowvec>(p-1);
    for(int j=0; j<m; j++) {
      if(Y_compl[j]>=t[i]){
        double XZbj = beta[0]*X_compl[j]+as_scalar(Z_compl.row(j)*beta.subvec(1,p-1));
        tempX += ipw_compl[j]*X_compl[j]*exp(XZbj);
        tempZ += ipw_compl[j]*Z_compl.row(j)*exp(XZbj);
      }
    }
    res(i,0) = tempX/sipw;
    res.submat(i,1,i,p-1) = tempZ/sipw;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube cpp_cox_joint_s2_org(arma::colvec beta, arma::colvec t, arma::mat data_compl, arma::colvec ipw_compl) {
  int p = beta.n_elem;
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec Y_compl = data_compl.col(1);
  arma::colvec X_compl = data_compl.col(3);
  arma::mat Z_compl = data_compl.submat(0,5,m-1,p+3);
  
  double sipw = sum(ipw_compl);
  
  arma::cube res(n,p,p);
  for(int i=0; i<n; i++) {
    double tempXX = 0;
    arma::rowvec tempXZ = zeros<arma::rowvec>(p-1);
    arma::mat tempZZ = zeros(p-1,p-1);
    for(int j=0; j<m; j++) {
      if(Y_compl[j]>=t[i]){
        double XZbj = beta[0]*X_compl[j]+as_scalar(Z_compl.row(j)*beta.subvec(1,p-1));
        tempXX += ipw_compl[j]*pow(X_compl[j],2)*exp(XZbj);
        tempXZ += ipw_compl[j]*X_compl[j]*Z_compl.row(j)*exp(XZbj);
        tempZZ += ipw_compl[j]*trans(Z_compl.row(j))*Z_compl.row(j)*exp(XZbj);
      }
    }
    arma::mat temp = zeros(p,p);
    temp(0,0) = tempXX;
    temp.submat(0,1,0,p-1) = tempXZ;
    temp.submat(1,0,p-1,0) = trans(tempXZ);
    temp.submat(1,1,p-1,p-1) = tempZZ;
    res.row(i) = temp/sipw;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_joint_s0_work(arma::colvec beta, arma::colvec t, arma::mat data_all) {
  int p = beta.n_elem;
  int n = t.n_elem;
  int m = data_all.n_rows;
  
  arma::colvec Y_all = data_all.col(1);
  arma::colvec X_star_all = data_all.col(4);
  arma::mat Z_all = data_all.submat(0,5,m-1,p+3);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp = 0;
    for(int j=0; j<m; j++) {
      if(Y_all[j]>=t[i]){
        temp += exp(beta[0]*X_star_all[j]+as_scalar(Z_all.row(j)*beta.subvec(1,p-1)));
      }
    }
    res[i] = temp/m;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_s1_work(arma::colvec beta, arma::colvec t, arma::mat data_all) {
  int p = beta.n_elem;
  int n = t.n_elem;
  int m = data_all.n_rows;
  
  arma::colvec Y_all = data_all.col(1);
  arma::colvec X_star_all = data_all.col(4);
  arma::mat Z_all = data_all.submat(0,5,m-1,p+3);
  
  arma::mat res = zeros(n,p);
  for(int i=0; i<n; i++) {
    double tempX = 0;
    arma::rowvec tempZ = zeros<arma::rowvec>(p-1);
    for(int j=0; j<m; j++) {
      if(Y_all[j]>=t[i]){
        double XZbj = beta[0]*X_star_all[j]+as_scalar(Z_all.row(j)*beta.subvec(1,p-1));
        tempX += X_star_all[j]*exp(XZbj);
        tempZ += Z_all.row(j)*exp(XZbj);
      }
    }
    res(i,0) = tempX/m;
    res.submat(i,1,i,p-1) = tempZ/m;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube cpp_cox_joint_s2_work(arma::colvec beta, arma::colvec t, arma::mat data_all) {
  int p = beta.n_elem;
  int n = t.n_elem;
  int m = data_all.n_rows;
  
  arma::colvec Y_all = data_all.col(1);
  arma::colvec X_star_all = data_all.col(4);
  arma::mat Z_all = data_all.submat(0,5,m-1,p+3);
  
  arma::cube res(n,p,p);
  for(int i=0; i<n; i++) {
    double tempXX = 0;
    arma::rowvec tempXZ = zeros<arma::rowvec>(p-1);
    arma::mat tempZZ = zeros(p-1,p-1);
    for(int j=0; j<m; j++) {
      if(Y_all[j]>=t[i]){
        double XZbj = beta[0]*X_star_all[j]+as_scalar(Z_all.row(j)*beta.subvec(1,p-1));
        tempXX += pow(X_star_all[j],2)*exp(XZbj);
        tempXZ += X_star_all[j]*Z_all.row(j)*exp(XZbj);
        tempZZ += trans(Z_all.row(j))*Z_all.row(j)*exp(XZbj);
      }
    }
    arma::mat temp = zeros(p,p);
    temp(0,0) = tempXX;
    temp.submat(0,1,0,p-1) = tempXZ;
    temp.submat(1,0,p-1,0) = trans(tempXZ);
    temp.submat(1,1,p-1,p-1) = tempZZ;
    res.row(i) = temp/m;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_joint_score_work_modelX(arma::colvec beta, arma::mat data, arma::mat data_all, arma::mat data_compl, arma::colvec predX, double seX, arma::colvec grid, arma::colvec weight, double pi, arma::colvec est_compl, arma::colvec ipw_compl, arma::colvec ipw) {
  int p = beta.n_elem;
  int n = data.n_rows;
  int n_all = data_all.n_rows;
  int q = grid.n_elem;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::mat Z = data.submat(0,5,n-1,p+3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_joint_s0_org(beta,Y,data_compl,ipw_compl);
  arma::mat s1 = cpp_cox_joint_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_joint_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::mat s1_all = cpp_cox_joint_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  double temp1X = 0;
  arma::rowvec temp1Z = zeros<arma::rowvec>(p-1);
  for(int i=0; i<n; i++) {
    double Lambda_Yi = cpp_cox_joint_est_Lambda(Y[i],data_compl,est_compl,ipw_compl);
    double Zbi = as_scalar(Z.row(i)*beta.subvec(1,p-1));
    
    double b = 0;
    double score = 0;
    double inf = 0;
    double diff_b = 10;
    for(int j=0; j<100; j++) {
      double old_b = b;
      score = D[i]*beta[0]-Lambda_Yi*exp(beta[0]*b+Zbi)*beta[0]+cpp_cox_joint_d_log_normal_pdf(b,predX[i],seX);
      inf = -Lambda_Yi*exp(beta[0]*b+Zbi)*pow(beta[0],2)+cpp_cox_joint_dd_log_normal_pdf(seX);
      b = b-score/inf;
      diff_b = abs(b-old_b);
      if(diff_b<1e-3) break;
    }
    arma::colvec gridb = grid/pow(-inf,1/2)+b;
    
    arma::colvec log_fb = zeros(q);
    for(int j=0; j<q; j++) {
      double log_fY_XZj = D[i]*(beta[0]*gridb[j]+Zbi)-Lambda_Yi*exp(beta[0]*gridb[j]+Zbi);
      log_fb[j] = log_fY_XZj+cpp_cox_joint_log_normal_pdf(gridb[j],predX[i],seX,pi);
    }
    double max_log_fb = max(log_fb);
    arma::colvec weight1 = zeros(q);
    double sum_weight1 = 0;
    for(int j=0; j<q; j++) {
      weight1[j] = exp(log_fb[j]-max_log_fb+pow(grid[j],2)+log(weight[j]));
      sum_weight1 += weight1[j];
    }
    arma::colvec weightb = weight1/sum_weight1;
    
    arma::colvec temp2X = D[i]*(gridb-s1(i,0)/s0[i]);
    arma::rowvec temp2Z = D[i]*(Z.row(i)-s1.submat(i,1,i,p-1)/s0[i]);
    arma::colvec temp3X = zeros(q);
    arma::mat temp3Z = zeros(q,p-1);
    for(int j=0; j<q; j++) {
      double temp4X = 0;
      arma::rowvec temp4Z = zeros<arma::rowvec>(p-1);
      for(int l=0; l<n_all; l++) {
        if(Y_all[l]<=Y[i]) {
          temp4X += D_all[l]*(gridb[j]-s1_all(l,0)/s0_all[l])/s0_all[l];
          temp4Z += D_all[l]*(Z.row(i)-s1_all.submat(l,1,l,p-1)/s0_all[l])/s0_all[l];
        }
      }
      temp3X[j] = exp(beta[0]*gridb[j]+Zbi)*temp4X/n_all;
      temp3Z.row(j) = exp(beta[0]*gridb[j]+Zbi)*temp4Z/n_all;
    }
    double temp5X = 0;
    arma::rowvec temp5Z = zeros<arma::rowvec>(p-1);
    for(int j=0; j<q; j++) {
      temp5X += (temp2X[j]-temp3X[j])*weightb[j];
      temp5Z += (temp2Z-temp3Z.row(j))*weightb[j];
    }
    temp1X += ipw[i]*temp5X;
    temp1Z += ipw[i]*temp5Z;
  }

  arma::colvec res = zeros(p);
  res[0] = temp1X;
  res.subvec(1,p-1) = trans(temp1Z);
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_score_work_modelX_vec(arma::colvec beta, arma::mat data, arma::mat data_all, arma::mat data_compl, arma::colvec predX, double seX, arma::colvec grid, arma::colvec weight, double pi, arma::colvec est_compl, arma::colvec ipw_compl, arma::colvec ipw) {
  int p = beta.n_elem;
  int n = data.n_rows;
  int n_all = data_all.n_rows;
  int q = grid.n_elem;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::mat Z = data.submat(0,5,n-1,p+3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_joint_s0_org(beta,Y,data_compl,ipw_compl);
  arma::mat s1 = cpp_cox_joint_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_joint_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::mat s1_all = cpp_cox_joint_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  arma::mat res = zeros(n,p);
  for(int i=0; i<n; i++) {
    double Lambda_Yi = cpp_cox_joint_est_Lambda(Y[i],data_compl,est_compl,ipw_compl);
    double Zbi = as_scalar(Z.row(i)*beta.subvec(1,p-1));
    
    double b = 0;
    double score = 0;
    double inf = 0;
    double diff_b = 10;
    for(int j=0; j<100; j++) {
      double old_b = b;
      score = D[i]*beta[0]-Lambda_Yi*exp(beta[0]*b+Zbi)*beta[0]+cpp_cox_joint_d_log_normal_pdf(b,predX[i],seX);
      inf = -Lambda_Yi*exp(beta[0]*b+Zbi)*pow(beta[0],2)+cpp_cox_joint_dd_log_normal_pdf(seX);
      b = b-score/inf;
      diff_b = abs(b-old_b);
      if(diff_b<1e-3) break;
    }
    arma::colvec gridb = grid/pow(-inf,1/2)+b;
    
    arma::colvec log_fb = zeros(q);
    for(int j=0; j<q; j++) {
      double log_fY_XZj = D[i]*(beta[0]*gridb[j]+Zbi)-Lambda_Yi*exp(beta[0]*gridb[j]+Zbi);
      log_fb[j] = log_fY_XZj+cpp_cox_joint_log_normal_pdf(gridb[j],predX[i],seX,pi);
    }
    double max_log_fb = max(log_fb);
    arma::colvec weight1 = zeros(q);
    double sum_weight1 = 0;
    for(int j=0; j<q; j++) {
      weight1[j] = exp(log_fb[j]-max_log_fb+pow(grid[j],2)+log(weight[j]));
      sum_weight1 += weight1[j];
    }
    arma::colvec weightb = weight1/sum_weight1;
    
    arma::colvec temp2X = D[i]*(gridb-s1(i,0)/s0[i]);
    arma::rowvec temp2Z = D[i]*(Z.row(i)-s1.submat(i,1,i,p-1)/s0[i]);
    arma::colvec temp3X = zeros(q);
    arma::mat temp3Z = zeros(q,p-1);
    for(int j=0; j<q; j++) {
      double temp4X = 0;
      arma::rowvec temp4Z = zeros<arma::rowvec>(p-1);
      for(int l=0; l<n_all; l++) {
        if(Y_all[l]<=Y[i]) {
          temp4X += D_all[l]*(gridb[j]-s1_all(l,0)/s0_all[l])/s0_all[l];
          temp4Z += D_all[l]*(Z.row(i)-s1_all.submat(l,1,l,p-1)/s0_all[l])/s0_all[l];
        }
      }
      temp3X[j] = exp(beta[0]*gridb[j]+Zbi)*temp4X/n_all;
      temp3Z.row(j) = exp(beta[0]*gridb[j]+Zbi)*temp4Z/n_all;
    }
    double temp5X = 0;
    arma::rowvec temp5Z = zeros<arma::rowvec>(p-1);
    for(int j=0; j<q; j++) {
      temp5X += (temp2X[j]-temp3X[j])*weightb[j];
      temp5Z += (temp2Z-temp3Z.row(j))*weightb[j];
    }
    res(i,0) = ipw[i]*temp5X;
    res.submat(i,1,i,p-1) = ipw[i]*temp5Z;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_joint_score_work_trueX(arma::colvec beta, arma::mat data, arma::mat data_all, arma::mat data_compl, arma::colvec grid, arma::colvec weight, double sde, double theta, double sdxz, double pi, arma::colvec ipw_compl, arma::colvec ipw) {
  int p = beta.n_elem;
  int n = data.n_rows;
  int n_all = data_all.n_rows;
  int q = grid.n_elem;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::mat Z = data.submat(0,5,n-1,p+3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_joint_s0_org(beta,Y,data_compl,ipw_compl);
  arma::mat s1 = cpp_cox_joint_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_joint_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::mat s1_all = cpp_cox_joint_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  double temp1X = 0;
  arma::rowvec temp1Z = zeros<arma::rowvec>(p-1);
  for(int i=0; i<n; i++) {
    double Lambda_Yi_true = cpp_cox_joint_Lambda(Y[i]);
    double Zbi = as_scalar(Z.row(i)*beta.subvec(1,p-1));
    
    double b = 0;
    double score = 0;
    double inf = 0;
    double diff_b = 10;
    for(int j=0; j<100; j++) {
      double old_b = b;
      score = D[i]*beta[0]-Lambda_Yi_true*exp(beta[0]*b+Zbi)*beta[0]+cpp_cox_joint_dm_log_normal_pdf(X_star[i],b,sde)+cpp_cox_joint_d_log_normal_pdf(b,theta*sum(Z.row(i)),sdxz);
      inf = -Lambda_Yi_true*exp(beta[0]*b+Zbi)*pow(beta[0],2)+cpp_cox_joint_ddm_log_normal_pdf(sde)+cpp_cox_joint_dd_log_normal_pdf(sdxz);
      b = b-score/inf;
      diff_b = abs(b-old_b);
      if(diff_b<1e-3) break;
    }
    arma::colvec gridb = grid/pow(-inf,1/2)+b;
    
    arma::colvec log_fb = zeros(q);
    for(int j=0; j<q; j++) {
      double log_fY_XZj = D[i]*(beta[0]*gridb[j]+Zbi)-Lambda_Yi_true*exp(beta[0]*gridb[j]+Zbi);
      log_fb[j] = log_fY_XZj+cpp_cox_joint_log_normal_pdf(X_star[i],gridb[j],sde,pi)+cpp_cox_joint_log_normal_pdf(gridb[j],theta*sum(Z.row(i)),sdxz,pi);
    }
    double max_log_fb = max(log_fb);
    arma::colvec weight1 = zeros(q);
    double sum_weight1 = 0;
    for(int j=0; j<q; j++) {
      weight1[j] = exp(log_fb[j]-max_log_fb+pow(grid[j],2)+log(weight[j]));
      sum_weight1 += weight1[j];
    }
    arma::colvec weightb = weight1/sum_weight1;
    
    arma::colvec temp2X = D[i]*(gridb-s1(i,0)/s0[i]);
    arma::rowvec temp2Z = D[i]*(Z.row(i)-s1.submat(i,1,i,p-1)/s0[i]);
    arma::colvec temp3X = zeros(q);
    arma::mat temp3Z = zeros(q,p-1);
    for(int j=0; j<q; j++) {
      double temp4X = 0;
      arma::rowvec temp4Z = zeros<arma::rowvec>(p-1);
      for(int l=0; l<n_all; l++) {
        if(Y_all[l]<=Y[i]) {
          temp4X += D_all[l]*(gridb[j]-s1_all(l,0)/s0_all[l])/s0_all[l];
          temp4Z += D_all[l]*(Z.row(i)-s1_all.submat(l,1,l,p-1)/s0_all[l])/s0_all[l];
        }
      }
      temp3X[j] = exp(beta[0]*gridb[j]+Zbi)*temp4X/n_all;
      temp3Z.row(j) = exp(beta[0]*gridb[j]+Zbi)*temp4Z/n_all;
    }
    double temp5X = 0;
    arma::rowvec temp5Z = zeros<arma::rowvec>(p-1);
    for(int j=0; j<q; j++) {
      temp5X += (temp2X[j]-temp3X[j])*weightb[j];
      temp5Z += (temp2Z-temp3Z.row(j))*weightb[j];
    }
    temp1X += ipw[i]*temp5X;
    temp1Z += ipw[i]*temp5Z;
  }
  
  arma::colvec res = zeros(p);
  res[0] = temp1X;
  res.subvec(1,p-1) = trans(temp1Z);
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_score_work_trueX_vec(arma::colvec beta, arma::mat data, arma::mat data_all, arma::mat data_compl, arma::colvec grid, arma::colvec weight, double sde, double theta, double sdxz, double pi, arma::colvec ipw_compl, arma::colvec ipw) {
  int p = beta.n_elem;
  int n = data.n_rows;
  int n_all = data_all.n_rows;
  int q = grid.n_elem;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::mat Z = data.submat(0,5,n-1,p+3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_joint_s0_org(beta,Y,data_compl,ipw_compl);
  arma::mat s1 = cpp_cox_joint_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_joint_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::mat s1_all = cpp_cox_joint_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  arma::mat res = zeros(n,p);
  for(int i=0; i<n; i++) {
    double Lambda_Yi_true = cpp_cox_joint_Lambda(Y[i]);
    double Zbi = as_scalar(Z.row(i)*beta.subvec(1,p-1));
    
    double b = 0;
    double score = 0;
    double inf = 0;
    double diff_b = 10;
    for(int j=0; j<100; j++) {
      double old_b = b;
      score = D[i]*beta[0]-Lambda_Yi_true*exp(beta[0]*b+Zbi)*beta[0]+cpp_cox_joint_dm_log_normal_pdf(X_star[i],b,sde)+cpp_cox_joint_d_log_normal_pdf(b,theta*sum(Z.row(i)),sdxz);
      inf = -Lambda_Yi_true*exp(beta[0]*b+Zbi)*pow(beta[0],2)+cpp_cox_joint_ddm_log_normal_pdf(sde)+cpp_cox_joint_dd_log_normal_pdf(sdxz);
      b = b-score/inf;
      diff_b = abs(b-old_b);
      if(diff_b<1e-3) break;
    }
    arma::colvec gridb = grid/pow(-inf,1/2)+b;
    
    arma::colvec log_fb = zeros(q);
    for(int j=0; j<q; j++) {
      double log_fY_XZj = D[i]*(beta[0]*gridb[j]+Zbi)-Lambda_Yi_true*exp(beta[0]*gridb[j]+Zbi);
      log_fb[j] = log_fY_XZj+cpp_cox_joint_log_normal_pdf(X_star[i],gridb[j],sde,pi)+cpp_cox_joint_log_normal_pdf(gridb[j],theta*sum(Z.row(i)),sdxz,pi);
    }
    double max_log_fb = max(log_fb);
    arma::colvec weight1 = zeros(q);
    double sum_weight1 = 0;
    for(int j=0; j<q; j++) {
      weight1[j] = exp(log_fb[j]-max_log_fb+pow(grid[j],2)+log(weight[j]));
      sum_weight1 += weight1[j];
    }
    arma::colvec weightb = weight1/sum_weight1;
    
    arma::colvec temp2X = D[i]*(gridb-s1(i,0)/s0[i]);
    arma::rowvec temp2Z = D[i]*(Z.row(i)-s1.submat(i,1,i,p-1)/s0[i]);
    arma::colvec temp3X = zeros(q);
    arma::mat temp3Z = zeros(q,p-1);
    for(int j=0; j<q; j++) {
      double temp4X = 0;
      arma::rowvec temp4Z = zeros<arma::rowvec>(p-1);
      for(int l=0; l<n_all; l++) {
        if(Y_all[l]<=Y[i]) {
          temp4X += D_all[l]*(gridb[j]-s1_all(l,0)/s0_all[l])/s0_all[l];
          temp4Z += D_all[l]*(Z.row(i)-s1_all.submat(l,1,l,p-1)/s0_all[l])/s0_all[l];
        }
      }
      temp3X[j] = exp(beta[0]*gridb[j]+Zbi)*temp4X/n_all;
      temp3Z.row(j) = exp(beta[0]*gridb[j]+Zbi)*temp4Z/n_all;
    }
    double temp5X = 0;
    arma::rowvec temp5Z = zeros<arma::rowvec>(p-1);
    for(int j=0; j<q; j++) {
      temp5X += (temp2X[j]-temp3X[j])*weightb[j];
      temp5Z += (temp2Z-temp3Z.row(j))*weightb[j];
    }
    res(i,0) = ipw[i]*temp5X;
    res.submat(i,1,i,p-1) = ipw[i]*temp5Z;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_joint_score_org(arma::colvec beta, arma::mat data, arma::mat data_all, arma::mat data_compl, arma::colvec ipw_compl, arma::colvec ipw) {
  int p = beta.n_elem;
  int n = data.n_rows;
  int n_all = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X = data.col(3);
  arma::mat Z = data.submat(0,5,n-1,p+3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_joint_s0_org(beta,Y,data_compl,ipw_compl);
  arma::mat s1 = cpp_cox_joint_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_joint_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::mat s1_all = cpp_cox_joint_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  double temp1X = 0;
  arma::rowvec temp1Z = zeros<arma::rowvec>(p-1);
  for(int i=0; i<n; i++) {
    double XZbi = beta[0]*X[i]+as_scalar(Z.row(i)*beta.subvec(1,p-1));
    double temp2X = D[i]*(X[i]-s1(i,0)/s0[i]);
    arma::rowvec temp2Z = D[i]*(Z.row(i)-s1.submat(i,1,i,p-1)/s0[i]);
    double temp3X = 0;
    arma::rowvec temp3Z = zeros<arma::rowvec>(p-1);
    for(int l=0; l<n_all; l++) {
      if(Y_all[l]<=Y[i]) {
        temp3X += D_all[l]*exp(XZbi)*(X[i]-s1_all(l,0)/s0_all[l])/s0_all[l];
        temp3Z += D_all[l]*exp(XZbi)*(Z.row(i)-s1_all.submat(l,1,l,p-1)/s0_all[l])/s0_all[l];
      }
    }
    temp1X += ipw[i]*(temp2X-temp3X/n_all);
    temp1Z += ipw[i]*(temp2Z-temp3Z/n_all);
  }
  
  arma::colvec res = zeros(p);
  res[0] = temp1X;
  res.subvec(1,p-1) = trans(temp1Z);
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_score_org_vec(arma::colvec beta, arma::mat data, arma::mat data_all, arma::mat data_compl, arma::colvec ipw_compl, arma::colvec ipw) {
  int p = beta.n_elem;
  int n = data.n_rows;
  int n_all = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X = data.col(3);
  arma::mat Z = data.submat(0,5,n-1,p+3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_joint_s0_org(beta,Y,data_compl,ipw_compl);
  arma::mat s1 = cpp_cox_joint_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_joint_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::mat s1_all = cpp_cox_joint_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  arma::mat res = zeros(n,p);
  for(int i=0; i<n; i++) {
    double XZbi = beta[0]*X[i]+as_scalar(Z.row(i)*beta.subvec(1,p-1));
    double temp2X = D[i]*(X[i]-s1(i,0)/s0[i]);
    arma::rowvec temp2Z = D[i]*(Z.row(i)-s1.submat(i,1,i,p-1)/s0[i]);
    double temp3X = 0;
    arma::rowvec temp3Z = zeros<arma::rowvec>(p-1);
    for(int l=0; l<n_all; l++) {
      if(Y_all[l]<=Y[i]) {
        temp3X += D_all[l]*exp(XZbi)*(X[i]-s1_all(l,0)/s0_all[l])/s0_all[l];
        temp3Z += D_all[l]*exp(XZbi)*(Z.row(i)-s1_all.submat(l,1,l,p-1)/s0_all[l])/s0_all[l];
      }
    }
    res(i,0) = ipw[i]*(temp2X-temp3X/n_all);
    res.submat(i,1,i,p-1) = ipw[i]*(temp2Z-temp3Z/n_all);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_info_org(arma::colvec beta, arma::mat data_all, arma::mat data_compl, arma::colvec ipw_compl) {
  int p = beta.n_elem;
  int n_all = data_all.n_rows;
  
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0_all = cpp_cox_joint_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::mat s1_all = cpp_cox_joint_s1_org(beta,Y_all,data_compl,ipw_compl);
  arma::cube s2_all = cpp_cox_joint_s2_org(beta,Y_all,data_compl,ipw_compl);
  
  arma::mat res = zeros(p,p);
  for(int i=0; i<n_all; i++) {
    arma::mat temp1 = s2_all.row(i);
    arma::mat temp2 = trans(s1_all.row(i))*s1_all.row(i);
    res += D_all[i]*(temp1/s0_all[i]-temp2/pow(s0_all[i],2));
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_joint_score_work_coxph(arma::colvec beta, arma::mat data, arma::mat data_all, arma::colvec ipw) {
  int p = beta.n_elem;
  int n = data.n_rows;
  int n_all = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::mat Z = data.submat(0,5,n-1,p+3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_joint_s0_work(beta,Y,data_all);
  arma::mat s1 = cpp_cox_joint_s1_work(beta,Y,data_all);
  
  arma::colvec s0_all = cpp_cox_joint_s0_work(beta,Y_all,data_all);
  arma::mat s1_all = cpp_cox_joint_s1_work(beta,Y_all,data_all);
  
  double temp1X = 0;
  arma::rowvec temp1Z = zeros<arma::rowvec>(p-1);
  for(int i=0; i<n; i++) {
    double XZbi = beta[0]*X_star[i]+as_scalar(Z.row(i)*beta.subvec(1,p-1));
    double temp2X = D[i]*(X_star[i]-s1(i,0)/s0[i]);
    arma::rowvec temp2Z = D[i]*(Z.row(i)-s1.submat(i,1,i,p-1)/s0[i]);
    double temp3X = 0;
    arma::rowvec temp3Z = zeros<arma::rowvec>(p-1);
    for(int l=0; l<n_all; l++) {
      if(Y_all[l]<=Y[i]) {
        temp3X += D_all[l]*exp(XZbi)*(X_star[i]-s1_all(l,0)/s0_all[l])/s0_all[l];
        temp3Z += D_all[l]*exp(XZbi)*(Z.row(i)-s1_all.submat(l,1,l,p-1)/s0_all[l])/s0_all[l];
      }
    }
    temp1X += ipw[i]*(temp2X-temp3X/n_all);
    temp1Z += ipw[i]*(temp2Z-temp3Z/n_all);
  }
  
  arma::colvec res = zeros(p);
  res[0] = temp1X;
  res.subvec(1,p-1) = trans(temp1Z);
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_score_work_coxph_vec(arma::colvec beta, arma::mat data, arma::mat data_all, arma::colvec ipw) {
  int p = beta.n_elem;
  int n = data.n_rows;
  int n_all = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::mat Z = data.submat(0,5,n-1,p+3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_joint_s0_work(beta,Y,data_all);
  arma::mat s1 = cpp_cox_joint_s1_work(beta,Y,data_all);
  
  arma::colvec s0_all = cpp_cox_joint_s0_work(beta,Y_all,data_all);
  arma::mat s1_all = cpp_cox_joint_s1_work(beta,Y_all,data_all);
  
  arma::mat res = zeros(n,p);
  for(int i=0; i<n; i++) {
    double XZbi = beta[0]*X_star[i]+as_scalar(Z.row(i)*beta.subvec(1,p-1));
    double temp2X = D[i]*(X_star[i]-s1(i,0)/s0[i]);
    arma::rowvec temp2Z = D[i]*(Z.row(i)-s1.submat(i,1,i,p-1)/s0[i]);
    double temp3X = 0;
    arma::rowvec temp3Z = zeros<arma::rowvec>(p-1);
    for(int l=0; l<n_all; l++) {
      if(Y_all[l]<=Y[i]) {
        temp3X += D_all[l]*exp(XZbi)*(X_star[i]-s1_all(l,0)/s0_all[l])/s0_all[l];
        temp3Z += D_all[l]*exp(XZbi)*(Z.row(i)-s1_all.submat(l,1,l,p-1)/s0_all[l])/s0_all[l];
      }
    }
    res(i,0) = ipw[i]*(temp2X-temp3X/n_all);
    res.submat(i,1,i,p-1) = ipw[i]*(temp2Z-temp3Z/n_all);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_joint_info_work_coxph(arma::colvec beta, arma::mat data_all) {
  int p = beta.n_elem;
  int n_all = data_all.n_rows;
  
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0_all = cpp_cox_joint_s0_work(beta,Y_all,data_all);
  arma::mat s1_all = cpp_cox_joint_s1_work(beta,Y_all,data_all);
  arma::cube s2_all = cpp_cox_joint_s2_work(beta,Y_all,data_all);
  
  arma::mat res = zeros(p,p);
  for(int i=0; i<n_all; i++) {
    arma::mat temp1 = s2_all.row(i);
    arma::mat temp2 = trans(s1_all.row(i))*s1_all.row(i);
    res += D_all[i]*(temp1/s0_all[i]-temp2/pow(s0_all[i],2));
  }
  
  return res;
}
