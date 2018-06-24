/*
 * Sample from a truncated mutivariate normal distribution
 * mu is K x 1; Sigma is K x K
 * returns a K x 1 vector
 */

inline
arma::vec
rmvnorm_trunc(const arma::vec& mu, const arma::mat& Sigma, const double a, const double b, bool pre_chol = false)
{
    //
    
    arma::vec ret = stats::rmvnorm(mu,Sigma,pre_chol);
    
	arma::uvec trunc_test_l = arma::find(ret < a);
	arma::uvec trunc_test_u = arma::find(ret > b);
	
	int check_1 = trunc_test_l.n_elem;
    int check_2 = trunc_test_u.n_elem;
    
    //
    
	if (check_1 + check_2 > 0) {
		while (check_1 + check_2 > 0) {
			ret = stats::rmvnorm(mu,Sigma,pre_chol);
			
			trunc_test_l = arma::find(ret < a);
			trunc_test_u = arma::find(ret > b);
			
			check_1 = trunc_test_l.n_elem;
			check_2 = trunc_test_u.n_elem;
		}
    }
    
    //
    
	return ret;
}
