#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat FARest(arma::mat& y,arma::vec& u, double u0,double p,double d,double bwp){

  arma::vec pd = {p,d};
  double lagrm = max(pd);
  int Tlength = y.n_rows;
  int k = y.n_cols;

  arma::mat mat_Y = y(arma::span(lagrm,Tlength-1),arma::span(0,k-1));

  arma::mat mat_X(Tlength-lagrm,k*p);
  for(int l = 0; l < p; l++){
    for(int j = 0; j < k; j++){
      mat_X(arma::span::all,k*l+j) = y(arma::span(lagrm-l-1,Tlength-l-2),j);
    }
  }

  arma::mat mat_U = u(arma::span(lagrm-d,Tlength-d-1));
  double bw = bwp * as_scalar(range(mat_U));
  arma::mat W = ((1/arma::datum::sqrt2pi)*exp(-0.5*((mat_U-u0)/bw)%((mat_U-u0)/bw)))/bw;
  arma::mat mat_W = diagmat(W);

  arma::mat mat_UX(Tlength-lagrm,k*p);
  for(int l = 0; l < k*p; l++){
    mat_UX(arma::span::all,l) = mat_X(arma::span::all,l)%(mat_U-u0);
  }

  arma::mat mat_Xcurl = join_rows(mat_X,mat_UX);
  arma::mat mat_XW = trans(mat_Xcurl) * mat_W;
  arma::mat mat_Q = mat_XW * mat_Xcurl;
  arma::mat mat_XWY = mat_XW * mat_Y;
  arma::mat mat_Qinv = inv(mat_Q);
  arma::mat fhat = mat_Qinv * mat_XWY;

  return fhat;
}

// [[Rcpp::export]]
arma::mat MXFARest(arma::vec SeriesNum,int Tlength,arma::mat& y,arma::vec& u,double u0,double p,double d,double bwp){

  arma::vec pd = {p,d};
  double lagrm = max(pd);
  int k = y.n_cols;
  int g_num = SeriesNum.n_elem;
  int SeriesNum_Total = sum(SeriesNum);
  arma::vec SeriesNum_cum = cumsum(SeriesNum);
  arma::vec SeriesNum_cumlag = SeriesNum_cum - SeriesNum;

  arma::cube fhat_comp(2*k*p,k,SeriesNum_Total);
  arma::cube fhat_comp_res(2*k*p,k,SeriesNum_Total);

  for(int i = 0; i < SeriesNum_Total; i++){
    arma::mat y_i = y(arma::span(i*Tlength,(i+1)*Tlength-1),arma::span::all);
    arma::vec u_i = u(arma::span(i*Tlength,(i+1)*Tlength-1));
    arma::mat fhat_comp_i = FARest(y_i,u_i,u0,p,d,bwp);
    fhat_comp.slice(i) = fhat_comp_i;
  }

  for(int s = 0; s < g_num; s++){

    arma::mat fhat_comp_m(2*k*p,k);
    for(int i = 0; i < 2*k*p; i++){
      for(int j = 0; j < k; j++){
        arma::vec sub = fhat_comp(arma::span(i),arma::span(j),arma::span(SeriesNum_cumlag(s),SeriesNum_cum(s) - 1));
        fhat_comp_m(i,j) = mean(sub);
      }
    }

    for(int i = SeriesNum_cumlag(s); i < SeriesNum_cum(s); i++){
      fhat_comp_res.slice(i) = fhat_comp.slice(i) - fhat_comp_m;
    }

  }

  arma::mat vhat_re(2*k*p,k);
  for(int i = 0; i < 2*k*p; i++){
    for(int j = 0; j < k; j++){
      arma::vec sub = fhat_comp_res(arma::span(i),arma::span(j),arma::span::all);
      vhat_re(i,j) = var(sub);
    }
  }

  arma::mat mat_Y(SeriesNum_Total*(Tlength-lagrm),k);
  for(int i = 0; i < SeriesNum_Total; i++){
    mat_Y(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),arma::span::all)
    = y(arma::span(i*Tlength+lagrm,(i+1)*Tlength-1),arma::span::all);
  }

  arma::mat mat_X(SeriesNum_Total*(Tlength-lagrm),k*p);
  for(int i = 0; i < SeriesNum_Total; i++){
    for(int l = 0; l < p; l++){
      for(int j = 0; j < k; j++){
        mat_X(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),k*l+j)
        = y(arma::span(i*Tlength+lagrm-(l+1),(i+1)*Tlength-(l+1)-1),j);
      }
    }
  }

  arma::mat mat_U(SeriesNum_Total*(Tlength-lagrm),1);
  for(int i = 0; i < SeriesNum_Total; i++){
    mat_U(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),0)
    = u(arma::span(i*Tlength+lagrm-d,(i+1)*Tlength-d-1));
  }

  arma::mat mat_UX(SeriesNum_Total*(Tlength-lagrm),k*p);
  for(int l = 0; l < k*p; l++){
    mat_UX(arma::span::all,l) = mat_X(arma::span::all,l)%(mat_U-u0);
  }

  arma::mat mat_Xcurl = join_rows(mat_X,mat_UX);

  arma::mat fhat(k,2*k*p*(SeriesNum_Total+g_num));

  for(int j = 0; j < k; j++){

    arma::mat fhat_j(2*k*p*(SeriesNum_Total+g_num),1);

    for(int s = 0; s < g_num; s++){

      arma::cube mat_A_all(2*k*p,2*k*p,SeriesNum(s));
      arma::cube mat_Omega_all(2*k*p,2*k*p,SeriesNum(s));
      arma::cube mat_R_all(2*k*p,1,SeriesNum(s));
      arma::cube mat_R2_all(2*k*p,1,SeriesNum(s));
      for(int i = SeriesNum_cumlag(s); i < SeriesNum_cum(s); i++){
        arma::mat mat_Xi = mat_Xcurl(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),arma::span::all);
        arma::mat mat_Ui = mat_U(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),0);
        double bw = bwp * as_scalar(range(mat_Ui));
        arma::mat Wi = ((1/arma::datum::sqrt2pi)*exp(-0.5*((mat_Ui-u0)/bw)%((mat_Ui-u0)/bw)))/bw;
        arma::mat mat_Wi = diagmat(Wi);
        double wgt = Wi.max();
        arma::mat sub = vhat_re(arma::span::all,j)/wgt;
        arma::mat mat_G = diagmat(sub);
        arma::mat mat_Yij = mat_Y(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),j);
        arma::mat mat_XWX = trans(mat_Xi) * mat_Wi * mat_Xi;
        arma::mat mat_XWY = trans(mat_Xi) * mat_Wi * mat_Yij;
        arma::mat mat_XWXG_inv = inv(trans(mat_Xi) * mat_Wi * mat_Xi + mat_G);
        arma::mat mat_Omega_i = mat_XWX - mat_XWX * mat_XWXG_inv * mat_XWX;
        arma::mat mat_R2_i = mat_XWY - mat_XWX * mat_XWXG_inv * mat_XWY;
        mat_A_all.slice(i-SeriesNum_cumlag(s)) =  mat_XWX;
        mat_Omega_all.slice(i-SeriesNum_cumlag(s)) = mat_Omega_i;
        mat_R_all.slice(i-SeriesNum_cumlag(s)) = mat_XWY;
        mat_R2_all.slice(i-SeriesNum_cumlag(s)) = mat_R2_i;
      }

      arma::mat mat_A(2*k*p,2*k*p);
      arma::mat mat_Omega(2*k*p,2*k*p);
      arma::mat mat_R(2*k*p,1);
      arma::mat mat_R2(2*k*p,1);
      for(int i = 0; i < 2*k*p; i++){
        for(int j = 0; j < 2*k*p; j++){
          arma::vec mat_A_sub = mat_A_all(arma::span(i),arma::span(j),arma::span::all);
          mat_A(i,j) = sum(mat_A_sub);
          arma::vec mat_Omega_sub = mat_Omega_all(arma::span(i),arma::span(j),arma::span::all);
          mat_Omega(i,j) = sum(mat_Omega_sub);
        }
        arma::vec mat_R_sub = mat_R_all(arma::span(i),arma::span(0),arma::span::all);
        mat_R(i,0) = sum(mat_R_sub);
        arma::vec mat_R2_sub = mat_R2_all(arma::span(i),arma::span(0),arma::span::all);
        mat_R2(i,0) = sum(mat_R2_sub);
      }

      arma::mat mat_A_inv = inv(mat_A);
      arma::mat mat_Omega_inv = inv(mat_Omega);
      arma::mat mat_Beta = mat_Omega_inv * mat_R2;

      arma::mat mat_Alpha(2*k*p*SeriesNum(s),1);
      for(int i = SeriesNum_cumlag(s); i < SeriesNum_cum(s); i++){
        arma::mat mat_Xi = mat_Xcurl(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),arma::span::all);
        arma::mat mat_Ui = mat_U(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),0);
        double bw = bwp * as_scalar(range(mat_Ui));
        arma::mat Wi = ((1/arma::datum::sqrt2pi)*exp(-0.5*((mat_Ui-u0)/bw)%((mat_Ui-u0)/bw)))/bw;
        arma::mat mat_Wi = diagmat(Wi);
        double wgt = Wi.max();
        arma::mat sub = wgt/vhat_re(arma::span::all,j);
        arma::mat mat_G = diagmat(sub);
        arma::mat mat_Yij = mat_Y(arma::span(i*(Tlength-lagrm),(i+1)*(Tlength-lagrm)-1),j);
        arma::mat mat_XWX = trans(mat_Xi) * mat_Wi * mat_Xi;
        arma::mat mat_XWY = trans(mat_Xi) * mat_Wi * mat_Yij;
        arma::mat mat_Phi_i = inv(mat_XWX + mat_G - mat_XWX * mat_A_inv * mat_XWX);
        arma::mat mat_R3_i = mat_XWY - mat_XWX * mat_A_inv * mat_R;
        arma::mat mat_Alpha_i = mat_Phi_i * mat_R3_i;
        mat_Alpha(arma::span((i-SeriesNum_cumlag(s))*(2*k*p),(i-SeriesNum_cumlag(s)+1)*(2*k*p)-1),0) = mat_Alpha_i;
      }

      fhat_j(arma::span(2*k*p*s,2*k*p*(s+1)-1),0) = mat_Beta;
      fhat_j(arma::span(2*k*p*(g_num + SeriesNum_cumlag(s)),2*k*p*(g_num + SeriesNum_cum(s))-1),0) = mat_Alpha;

    }

    fhat(j,arma::span::all) = trans(fhat_j);
  }

  return fhat;

}

