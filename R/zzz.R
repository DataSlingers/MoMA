.onAttach <- function(...) { # nocov start
    if(interactive()){
        MoMA::moma_set_logger_level_cpp(0)
        msg <- c("Thank you for using MoMA!",
                 "The current logging level is",
                 sQuote(paste0(moma_logger_level(), ".")),
                 "To change this, see ?moma_logging.")

        packageStartupMessage(paste(msg, collapse=" "))
    }
}                            # nocov end
