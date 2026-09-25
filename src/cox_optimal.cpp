#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_lambda(double t) {
  double res = 0.5*t;
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_Lambda(double t) {
  double res = 0.25*pow(t,2);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_log_kernel_fn(double t, double h, double pi) {
  double res = -pow(t/h,2)/2-log(2*pi)/2-log(h);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_log_normal_pdf(double t, double m, double s, double pi) {
  double res = -pow((t-m)/s,2)/2-log(2*pi)/2-log(s);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_est_Lambda(double t, arma::mat data_compl, double est_compl, arma::colvec ipw_compl) {
  int m = data_compl.n_rows;
  
  arma::colvec Y = data_compl.col(1);
  arma::colvec D = data_compl.col(2);
  arma::colvec X = data_compl.col(3);
  
  double res = 0;
  for(int i=0; i<m; i++) {
    if((Y[i]<=t)&&(D[i]==1)) {
      double temp = 0;
      for(int j=0; j<m; j++) {
        if(Y[j]>=Y[i]) {
          temp += ipw_compl[j]*exp(est_compl*X[j]);
        }
      }
      res += ipw_compl[i]/temp;
    }
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_est_lambda(arma::mat data_compl, double est_compl, arma::colvec ipw_compl) {
  int m = data_compl.n_rows;
  
  arma::colvec Y = data_compl.col(1);
  arma::colvec D = data_compl.col(2);
  arma::colvec X = data_compl.col(3);
  
  arma::colvec res = zeros(m);
  for(int i=0; i<m; i++) {
    if(D[i]==1) {
      double temp = 0;
      for(int j=0; j<m; j++) {
        if(Y[j]>=Y[i]) {
          temp += ipw_compl[j]*exp(est_compl*X[j]);
        }
      }
      res[i] = ipw_compl[i]/temp;
    }
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_s0_true(double beta, arma::colvec t, arma::mat data_compl, double est_compl, arma::colvec ipw_compl) {
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec X_compl = data_compl.col(3);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double Lambda_ti = cpp_cox_opt_est_Lambda(t[i],data_compl,est_compl,ipw_compl);
    double temp = 0;
    for(int j=0; j<m; j++) {
      temp += ipw_compl[j]*exp(beta*X_compl[j]-Lambda_ti*exp(beta*X_compl[j]));
    }
    res[i] = temp/sum(ipw_compl);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_s1_true(double beta, arma::colvec t, arma::mat data_compl, double est_compl, arma::colvec ipw_compl) {
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec X_compl = data_compl.col(3);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double Lambda_ti = cpp_cox_opt_est_Lambda(t[i],data_compl,est_compl,ipw_compl);
    double temp = 0;
    for(int j=0; j<m; j++) {
      temp += ipw_compl[j]*X_compl[j]*exp(beta*X_compl[j]-Lambda_ti*exp(beta*X_compl[j]));
    }
    res[i] = temp/sum(ipw_compl);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_s0_org(double beta, arma::colvec t, arma::mat data_compl, arma::colvec ipw_compl) {
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec Y_compl = data_compl.col(1);
  arma::colvec X_compl = data_compl.col(3);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp = 0;
    for(int j=0; j<m; j++) {
      if(Y_compl[j]>=t[i]) {
        temp += ipw_compl[j]*exp(beta*X_compl[j]);
      }
    }
    res[i] = temp/sum(ipw_compl);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_s1_org(double beta, arma::colvec t, arma::mat data_compl, arma::colvec ipw_compl) {
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec Y_compl = data_compl.col(1);
  arma::colvec X_compl = data_compl.col(3);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp = 0;
    for(int j=0; j<m; j++) {
      if(Y_compl[j]>=t[i]) {
        temp += ipw_compl[j]*X_compl[j]*exp(beta*X_compl[j]);
      }
    }
    res[i] = temp/sum(ipw_compl);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_s2_org(double beta, arma::colvec t, arma::mat data_compl, arma::colvec ipw_compl) {
  int n = t.n_elem;
  int m = data_compl.n_rows;
  
  arma::colvec Y_compl = data_compl.col(1);
  arma::colvec X_compl = data_compl.col(3);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp = 0;
    for(int j=0; j<m; j++) {
      if(Y_compl[j]>=t[i]) {
        temp += ipw_compl[j]*pow(X_compl[j],2)*exp(beta*X_compl[j]);
      }
    }
    res[i] = temp/sum(ipw_compl);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_s0_work(double beta, arma::colvec t, arma::mat data_all) {
  int n = t.n_elem;
  int m = data_all.n_rows;
  
  arma::colvec Y_all = data_all.col(1);
  arma::colvec X_star_all = data_all.col(4);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp = 0;
    for(int j=0; j<m; j++) {
      if(Y_all[j]>=t[i]) {
        temp += exp(beta*X_star_all[j]);
      }
    }
    res[i] = temp/m;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_s1_work(double beta, arma::colvec t, arma::mat data_all) {
  int n = t.n_elem;
  int m = data_all.n_rows;
  
  arma::colvec Y_all = data_all.col(1);
  arma::colvec X_star_all = data_all.col(4);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp = 0;
    for(int j=0; j<m; j++) {
      if(Y_all[j]>=t[i]) {
        temp += X_star_all[j]*exp(beta*X_star_all[j]);
      }
    }
    res[i] = temp/m;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_s2_work(double beta, arma::colvec t, arma::mat data_all) {
  int n = t.n_elem;
  int m = data_all.n_rows;
  
  arma::colvec Y_all = data_all.col(1);
  arma::colvec X_star_all = data_all.col(4);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp = 0;
    for(int j=0; j<m; j++) {
      if(Y_all[j]>=t[i]) {
        temp += pow(X_star_all[j],2)*exp(beta*X_star_all[j]);
      }
    }
    res[i] = temp/m;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_score_work_kernel(double beta, arma::mat data, arma::mat data_all, arma::mat data_compl, double h, double pi, double est_compl, arma::colvec ipw_compl, arma::colvec ipw) {
  int n = data.n_rows;
  int m = data_compl.n_rows;
  int n_all = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::colvec X_compl = data_compl.col(3);
  arma::colvec X_star_compl = data_compl.col(4);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_opt_s0_org(beta,Y,data_compl,ipw_compl);
  arma::colvec s1 = cpp_cox_opt_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_opt_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::colvec s1_all = cpp_cox_opt_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  double res = 0;
  for(int i=0; i<n; i++) {
    double Lambda_Yi = cpp_cox_opt_est_Lambda(Y[i],data_compl,est_compl,ipw_compl);
    arma::colvec temp1 = D[i]*(X_compl-s1[i]/s0[i]);
    arma::colvec temp2 = zeros(m);
    for(int j=0; j<m; j++) {
      double temp3 = 0;
      for(int l=0; l<n_all; l++) {
        if(Y_all[l]<=Y[i]) {
          temp3 += D_all[l]*exp(beta*X_compl[j])*(X_compl[j]-s1_all[l]/s0_all[l])/s0_all[l];
        }
      }
      temp2[j] = temp3/n_all;
    }
    
    arma::colvec log_fY_X = zeros(m);
    arma::colvec log_weight = zeros(m);
    for(int j=0; j<m; j++) {
      log_fY_X[j] = D[i]*beta*X_compl[j]-Lambda_Yi*exp(beta*X_compl[j]);
      log_weight[j] = cpp_cox_opt_log_kernel_fn(X_star[i]-X_star_compl[j], h, pi);
    }
    
    double num = 0;
    double den = 0;
    for(int j=0; j<m; j++) {
      double fY_Xj = exp(log_fY_X[j]-max(log_fY_X));
      double weightj = exp(log_weight[j]-max(log_weight));
      num += ipw_compl[j]*(temp1[j]-temp2[j])*fY_Xj*weightj;
      den += ipw_compl[j]*fY_Xj*weightj;
    }
    res += ipw[i]*num/den;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_score_work_kernel_vec(double beta, arma::mat data, arma::mat data_all, arma::mat data_compl, double h, double pi, double est_compl, arma::colvec ipw_compl, arma::colvec ipw) {
  int n = data.n_rows;
  int m = data_compl.n_rows;
  int n_all = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::colvec X_compl = data_compl.col(3);
  arma::colvec X_star_compl = data_compl.col(4);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_opt_s0_org(beta,Y,data_compl,ipw_compl);
  arma::colvec s1 = cpp_cox_opt_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_opt_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::colvec s1_all = cpp_cox_opt_s1_org(beta,Y_all,data_compl,ipw_compl);
  arma::colvec exp_beta_X = exp(beta*X_compl);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double Lambda_Yi = cpp_cox_opt_est_Lambda(Y[i],data_compl,est_compl,ipw_compl);
    arma::colvec temp1 = D[i]*(X_compl-s1[i]/s0[i]);
    double event_scale = 0;
    double event_mean_scale = 0;
    for(int l=0; l<n_all; l++) {
      if(Y_all[l]<=Y[i]) {
        double inv_s0 = 1/s0_all[l];
        event_scale += D_all[l]*inv_s0;
        event_mean_scale += D_all[l]*s1_all[l]*inv_s0*inv_s0;
      }
    }
    arma::colvec temp2 = exp_beta_X %
      (X_compl*event_scale-event_mean_scale)/n_all;
    
    arma::colvec log_fY_X = D[i]*beta*X_compl-Lambda_Yi*exp_beta_X;
    arma::colvec log_weight = zeros(m);
    for(int j=0; j<m; j++) {
      log_weight[j] = cpp_cox_opt_log_kernel_fn(X_star[i]-X_star_compl[j], h, pi);
    }
    double max_log_fY_X = max(log_fY_X);
    double max_log_weight = max(log_weight);
    
    double num = 0;
    double den = 0;
    for(int j=0; j<m; j++) {
      double fY_Xj = exp(log_fY_X[j]-max_log_fY_X);
      double weightj = exp(log_weight[j]-max_log_weight);
      num += ipw_compl[j]*(temp1[j]-temp2[j])*fY_Xj*weightj;
      den += ipw_compl[j]*fY_Xj*weightj;
    }
    res[i] = ipw[i]*num/den;
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_score_org(double beta, arma::mat data, arma::mat data_all, arma::mat data_compl, arma::colvec ipw_compl, arma::colvec ipw) {
  int n = data.n_rows;
  int m = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X = data.col(3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_opt_s0_org(beta,Y,data_compl,ipw_compl);
  arma::colvec s1 = cpp_cox_opt_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_opt_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::colvec s1_all = cpp_cox_opt_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  double res = 0;
  for(int i=0; i<n; i++) {
    double temp1 = D[i]*(X[i]-s1[i]/s0[i]);
    double temp2 = 0;
    for(int j=0; j<m; j++) {
      if(Y_all[j]<=Y[i]) {
        temp2 += D_all[j]*exp(beta*X[i])*(X[i]-s1_all[j]/s0_all[j])/s0_all[j];
      }
    }
    res += ipw[i]*(temp1-temp2/m);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_score_org_vec(double beta, arma::mat data, arma::mat data_all, arma::mat data_compl, arma::colvec ipw_compl, arma::colvec ipw) {
  int n = data.n_rows;
  int m = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X = data.col(3);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_opt_s0_org(beta,Y,data_compl,ipw_compl);
  arma::colvec s1 = cpp_cox_opt_s1_org(beta,Y,data_compl,ipw_compl);
  
  arma::colvec s0_all = cpp_cox_opt_s0_org(beta,Y_all,data_compl,ipw_compl);
  arma::colvec s1_all = cpp_cox_opt_s1_org(beta,Y_all,data_compl,ipw_compl);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp1 = D[i]*(X[i]-s1[i]/s0[i]);
    double temp2 = 0;
    for(int j=0; j<m; j++) {
      if(Y_all[j]<=Y[i]) {
        temp2 += D_all[j]*exp(beta*X[i])*(X[i]-s1_all[j]/s0_all[j])/s0_all[j];
      }
    }
    res[i] = ipw[i]*(temp1-temp2/m);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_info_org(double beta, arma::mat data_all, arma::mat data_compl, arma::colvec ipw_compl) {
  int m = data_all.n_rows;
  
  arma::colvec Y = data_all.col(1);
  arma::colvec D = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_opt_s0_org(beta,Y,data_compl,ipw_compl);
  arma::colvec s1 = cpp_cox_opt_s1_org(beta,Y,data_compl,ipw_compl);
  arma::colvec s2 = cpp_cox_opt_s2_org(beta,Y,data_compl,ipw_compl);
  
  double res = 0;
  for(int i=0; i<m; i++) {
    res += D[i]*(s2[i]/s0[i]-pow(s1[i]/s0[i],2));
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_score_work_coxph(double beta, arma::mat data, arma::mat data_all, arma::colvec ipw) {
  int n = data.n_rows;
  int m = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_opt_s0_work(beta,Y,data_all);
  arma::colvec s1 = cpp_cox_opt_s1_work(beta,Y,data_all);
  
  arma::colvec s0_all = cpp_cox_opt_s0_work(beta,Y_all,data_all);
  arma::colvec s1_all = cpp_cox_opt_s1_work(beta,Y_all,data_all);
  
  double res = 0;
  for(int i=0; i<n; i++) {
    double temp1 = D[i]*(X_star[i]-s1[i]/s0[i]);
    double temp2 = 0;
    for(int j=0; j<m; j++) {
      if(Y_all[j]<=Y[i]) {
        temp2 += D_all[j]*exp(beta*X_star[i])*(X_star[i]-s1_all[j]/s0_all[j])/s0_all[j];
      }
    }
    res += ipw[i]*(temp1-temp2/m);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec cpp_cox_opt_score_work_coxph_vec(double beta, arma::mat data, arma::mat data_all, arma::colvec ipw) {
  int n = data.n_rows;
  int m = data_all.n_rows;
  
  arma::colvec Y = data.col(1);
  arma::colvec D = data.col(2);
  arma::colvec X_star = data.col(4);
  arma::colvec Y_all = data_all.col(1);
  arma::colvec D_all = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_opt_s0_work(beta,Y,data_all);
  arma::colvec s1 = cpp_cox_opt_s1_work(beta,Y,data_all);
  
  arma::colvec s0_all = cpp_cox_opt_s0_work(beta,Y_all,data_all);
  arma::colvec s1_all = cpp_cox_opt_s1_work(beta,Y_all,data_all);
  
  arma::colvec res = zeros(n);
  for(int i=0; i<n; i++) {
    double temp1 = D[i]*(X_star[i]-s1[i]/s0[i]);
    double temp2 = 0;
    for(int j=0; j<m; j++) {
      if(Y_all[j]<=Y[i]) {
        temp2 += D_all[j]*exp(beta*X_star[i])*(X_star[i]-s1_all[j]/s0_all[j])/s0_all[j];
      }
    }
    res[i] = ipw[i]*(temp1-temp2/m);
  }
  
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double cpp_cox_opt_info_work_coxph(double beta, arma::mat data_all) {
  int m = data_all.n_rows;
  
  arma::colvec Y = data_all.col(1);
  arma::colvec D = data_all.col(2);
  
  arma::colvec s0 = cpp_cox_opt_s0_work(beta,Y,data_all);
  arma::colvec s1 = cpp_cox_opt_s1_work(beta,Y,data_all);
  arma::colvec s2 = cpp_cox_opt_s2_work(beta,Y,data_all);
  
  double res = 0;
  for(int i=0; i<m; i++) {
    res += D[i]*(s2[i]/s0[i]-pow(s1[i]/s0[i],2));
  }
  
  return res;
}
