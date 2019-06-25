
#ifndef MOMA_SOLVER_BICSEARCH_H
#define MOMA_SOLVER_BICSEARCH_H 1

#include "moma_base.h"
#include "moma_logging.h"
#include "moma_solver.h"

class BIC_searcher
{
  public:
    typedef double (PR_solver::*Criterion)(arma::vec y, const arma::vec &est);
    BIC_searcher(){};

    void bind(PR_solver *object, Criterion method);

    // current criterion
    double cur_criterion(arma::vec y, const arma::vec &est);

    ~BIC_searcher()
    {
        // No need to delete pr_solver
        MoMALogger::debug("Releasing a BIC_searcher object");
    }

    Rcpp::List search(const arma::vec &y,  // min_{u} || y - u || + ...penalty...
                      const arma::vec &u,  // start point
                      const arma::vec &alpha_u,
                      const arma::vec &lambda_u);

  private:
    PR_solver *pr_solver;
    Criterion cri;
};

#endif
