get_5Dlist_elem <- function(x, alpha_u_i, lambda_u_i, alpha_v_i, lambda_v_i, rank_i = 1) {
    if (!inherits(x, "MoMA_5D_list")) {
        moma_error(sQuote("x"), " should be a ", sQuote("MoMA_5D_list"), " object.")
    }
    n_alpha_u <- dim(x)[1]
    n_lambda_u <- dim(x)[2]
    n_alpha_v <- dim(x)[3]
    n_lambda_v <- dim(x)[4]
    n_rank <- dim(x)[5]

    # NOTE: R index starts from 1
    if (
        alpha_u_i <= 0 || alpha_u_i > n_alpha_u ||
            lambda_u_i <= 0 || lambda_u_i > n_lambda_u ||
            alpha_v_i <= 0 || alpha_v_i > n_alpha_v ||
            lambda_v_i <= 0 || lambda_v_i > n_lambda_v ||
            rank_i <= 0 || rank_i > n_rank
    ) {
        moma_error(
            "Invalid index (", alpha_u_i, ",", lambda_u_i,
            ",", alpha_v_i, ",", lambda_v_i, ",", rank_i, "), dim = ",
            dim(x)
        )
    }

    return(x[
        rank_i + n_rank * (
            lambda_v_i - 1 + n_lambda_v * (
                alpha_v_i - 1 + n_alpha_v * (
                    lambda_u_i - 1 + n_lambda_u * (
                        alpha_u_i - 1
                    )
                )
            )
        )
    ])
}
