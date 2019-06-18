#include "moma.h"

// [[Rcpp::export]]
void moma_set_logger_level_cpp(int level)
{
    auto logger_level = static_cast<MoMALoggerLevel>(level);
    MoMALogger::set_level(logger_level);
}

// [[Rcpp::export]]
int moma_get_logger_level_cpp()
{
    auto logger_level = static_cast<int>(MoMALogger::get_level());
    return logger_level;
}

// [[Rcpp::export]]
void moma_log_cpp(int level, Rcpp::StringVector x)
{
    auto msg_level  = static_cast<MoMALoggerLevel>(level);
    std::string msg = Rcpp::as<std::string>(x[0]);
    if (msg_level >= MoMALoggerLevel::ERRORS)
    {
        MoMALogger::error(msg);
    }
    else if (msg_level >= MoMALoggerLevel::WARNING)
    {
        MoMALogger::warning(msg);
    }
    else if (msg_level >= MoMALoggerLevel::MESSAGES)
    {
        MoMALogger::message(msg);
    }
    else if (msg_level >= MoMALoggerLevel::INFO)
    {
        MoMALogger::info(msg);
    }
    else if (msg_level >= MoMALoggerLevel::DEBUG)
    {
        MoMALogger::debug(msg);
    }
}
