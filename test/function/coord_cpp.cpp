// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; using namespace arma;


// [[Rcpp::export]]
arma::mat coord_cpp(const arma::colvec& tempSeas, const arma::rowvec& coeffTmp,
                const arma::colvec& CWD, const arma::rowvec& coeffCWD,
                const arma::colvec& PS, const arma::rowvec& coeffPS, 
                const arma::mat& WD_simu, const arma::mat& D_simu,
                const arma::rowvec& coeffWD, const arma::rowvec& coefflnD,
                const arma::rowvec& coefflnD2, const arma::rowvec& coeffE,
                const arma::rowvec& intercept) {
  
  
  arma::mat E_simu = tempSeas * coeffTmp + CWD * coeffCWD + PS * coeffPS;
  arma::mat l_WD_simu = log(WD_simu);
  arma::mat l_D_simu = log(D_simu);
  arma::mat l2_D_simu = square(l_D_simu);
  
  arma::mat AGB_simu = l_WD_simu.each_row() % coeffWD + l_D_simu.each_row() % coefflnD + l2_D_simu.each_row() % coefflnD2 - E_simu.each_row() % coeffE;
  AGB_simu = AGB_simu.each_row() + intercept;
  
  return AGB_simu;
}


/*** R
AGB_simu = coord_cpp(tempSeas = bioclimParams$tempSeas, coeffTmp = param_7$temp[selec], 
                     CWD = bioclimParams$CWD, coeffCWD = param_7$cwd[selec],
                     PS = bioclimParams$precSeas, coeffPS = param_7$prec[selec],
                     WD_simu = WD_simu, D_simu = D_simu, coeffWD =  param_7$logwsg[selec], coefflnD = param_7$logdbh[selec],
                     coefflnD2 = param_7$logdbh2[selec], coeffE = param_7$E[selec], intercept = param_7$intercept[selec])

AGB_simu1 = coord(n, coord12, WD_simu, D_simu, selec)
all( AGB_simu == AGB_simu1 )
*/
