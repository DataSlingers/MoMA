get_4Dlist_elem <- function(x, alpha_u_i, lambda_u_i, alpha_v_i, lambda_v_i){
    if(!inherits(x, "MoMA_4D_list")){
        moma_error(sQuote("x"), " should be a ", sQuote("MoMA_4D_list"), " object.")
    }
    n_alpha_u = dim(x)[1]
    n_lambda_u = dim(x)[2]
    n_alpha_v = dim(x)[3]
    n_lambda_v = dim(x)[4]

    # NOTE: R index starts from 1
    if(
        alpha_u_i  <= 0 || alpha_u_i > n_alpha_u ||
        lambda_u_i <= 0 || lambda_u_i > n_lambda_u ||
        alpha_v_i  <= 0 || alpha_v_i > n_alpha_v ||
        lambda_v_i <= 0 || lambda_v_i > n_lambda_v
    ){
        moma_error("Invalid index (",alpha_u_i, ",", lambda_u_i, 
                   ",",alpha_v_i, ",",lambda_v_i,"), dim = ",
                   dim(x))
    }
    return(x[n_lambda_u * n_alpha_v * n_lambda_v * (alpha_u_i-1) + 
                          n_alpha_v * n_lambda_v * (lambda_u_i-1) + 
                                      n_lambda_v * (alpha_v_i-1) +
                                                   lambda_v_i])
}
