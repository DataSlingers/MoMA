# include "moma.h"

class RcppFourDList{
    int n_au;
    int n_lu;
    int n_av;
    int n_lv;
    Rcpp::List my_list;

public:
    RcppFourDList(int au, int lu, int av, int lv):
        n_au(au), n_lu(lu), n_av(av), n_lv(lv), 
        my_list(n_au * n_av * n_lu * n_lv)
    {
            my_list.attr("dim") = Rcpp::NumericVector::create(
                n_au, n_av, n_lu, n_lv);
    };
 
    int insert(Rcpp::List object, int au_i, int lu_i, int av_i, int lv_i){
        // insert object in the au_i-th position along the au-axis
        // and so on
        if(
            au_i < 0 || au_i >= n_au ||
            lu_i < 0 || lu_i >= n_lu ||
            av_i < 0 || av_i >= n_av ||
            lv_i < 0 || lv_i >= n_lv
        ){
            MoMALogger::error("Invalid index is passed to RcppFourDList::insert. ")
                << "Dimension is (" << n_au << ", " 
                << n_lu << ", " << n_av << ", "
                << n_lv << "), (" 
                << au_i << ", " << lu_i << ", " << av_i << ", " << lv_i << ")."
                ;
        }
        my_list(n_lu * n_av * n_lv * au_i + 
                       n_av * n_lv * lu_i + 
                              n_lv * av_i +
                                     lv_i) = object;
        return 0;
    }

    Rcpp::List get_list(){
        return my_list;
    }
};
