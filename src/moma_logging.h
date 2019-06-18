#ifndef MOMA_LOGGING_H
#define MOMA_LOGGING_H 1

// Logging levels similar to Python -- see
// https://docs.python.org/3.6/library/logging.html#levels
// but changed ERROR ==> ERRORS to appease compiler and
// removed CRITICAL level.
//
// Also added a MESSAGES level for user-relevant info
// that is displayed by default (corresponding to message() in R.)
//
// This must be kept consistent with R/logging.R::MoMA_set_logger_level
enum class MoMALoggerLevel
{
    ERRORS   = 40,
    WARNING  = 30,
    MESSAGES = 20,
    INFO     = 10,
    DEBUG    = 0
};

// See http://gallery.rcpp.org/articles/quiet-stop-and-warning/
inline void warningNoCall(const std::string &s)
{
    Rf_warningcall(R_NilValue, s.c_str());
};

inline void NORET stopNoCall(const std::string &s)
{
    throw Rcpp::exception(s.c_str(), false);
};

// Logger structure loosely based on https://github.com/PennEcon/RcppLogger
class MoMALoggerMessage
{
  public:
    MoMALoggerMessage(const char *header,
                      MoMALoggerLevel msg_level,
                      MoMALoggerLevel logger_level,
                      std::ostream &logger_ostream)
    {
        this->msg_level    = msg_level;
        this->logger_level = logger_level;

        if (msg_level >= logger_level)
        {
            logger_ostream << header << " -- ";

            // Test for a reasonable compiler which supports std::put_time and
            // similar... If we have a recent _real_ GCC (i.e., GCC >= 5) or clang (>=
            // 4) then we should be good...
#if (__clang_major__ < 4) || ((__GNUC__ < 5) && !(__clang__))
#else
            auto time_now   = std::chrono::system_clock::now();
            auto time_now_t = std::chrono::system_clock::to_time_t(time_now);
            auto gmt_time   = gmtime(&time_now_t);

            logger_ostream << std::put_time(gmt_time, "%Y-%m-%d %H:%M:%S") << " -- ";
#endif
        }
    }

    ~MoMALoggerMessage()
    {
        if (msg_level >= logger_level)
        {
            logger_ostream << std::endl;
        }
    }

    template <typename T>
    MoMALoggerMessage &operator<<(const T &t)
    {
        if (msg_level >= logger_level)
        {
            logger_ostream << t;
        }
        return *this;
    }

  private:
    MoMALoggerLevel msg_level;
    MoMALoggerLevel logger_level;
    std::ostream &logger_ostream = Rcpp::Rcout;
};

// Special LoggerMessage for things we want to have
// handled via R's condition handling mechanisms
class RHandleMoMALoggerMessage
{
  public:
    RHandleMoMALoggerMessage(const char *header,
                             MoMALoggerLevel msg_level,
                             MoMALoggerLevel logger_level)
    {
        this->msg_level    = msg_level;
        this->logger_level = logger_level;
        this->log_msg      = new std::stringstream();
    }

    ~RHandleMoMALoggerMessage()
    {
        if (msg_level >= logger_level)
        {
            (*log_msg) << std::endl;
        }

        const std::string &log_msg_s = (*log_msg).str();
        delete log_msg;

        // Need an extra call to BEGIN_RCPP and VOID_END_RCPP
        // to handle exceptions thrown below
        // (they seem to interact badly with destructors otherwise)
        BEGIN_RCPP

        if (msg_level >= logger_level)
        {
            // Report to R (else: ignore)
            if (msg_level >= MoMALoggerLevel::ERRORS)
            {
                stopNoCall(log_msg_s.c_str());
            }
            else if (msg_level >= MoMALoggerLevel::WARNING)
            {
                warningNoCall(log_msg_s.c_str());
            }
            else if (msg_level >= MoMALoggerLevel::MESSAGES)
            {
                Rcpp::Function print_msg("message");
                print_msg(log_msg_s, Rcpp::Named("appendLF", false));
            }
        }

        VOID_END_RCPP
    }

    template <typename T>
    RHandleMoMALoggerMessage &operator<<(const T &t)
    {
        if (msg_level >= logger_level)
        {
            (*log_msg) << t;
        }
        return *this;
    }

    RHandleMoMALoggerMessage(RHandleMoMALoggerMessage &&) = default;
    RHandleMoMALoggerMessage &operator=(RHandleMoMALoggerMessage &&) = default;

  private:
    MoMALoggerLevel msg_level;
    MoMALoggerLevel logger_level;
    std::stringstream *log_msg;
};

// Singleton pattern loosely based on https://stackoverflow.com/a/1008289/967712
class MoMALogger
{
  public:
    static MoMALogger &getInstance()
    {
        // Guaranteed to be destroyed
        // Instantiated (below) on first use
        static MoMALogger instance;
        return instance;
    }

    // Delete these two methods but don't implement so we don't get default
    // implementations. This prevents "accidental" copies of the singleton
    // MoMALogger object.
    //
    // Note: Scott Meyers mentions in his Effective Modern
    //       C++ book, that deleted functions should generally
    //       be public as it results in better error messages
    //       due to the compilers behavior to check accessibility
    //       before deleted status
    MoMALogger(MoMALogger const &) = delete;      // Don't Implement
    void operator=(MoMALogger const &) = delete;  // Don't implement

    static RHandleMoMALoggerMessage error(const std::string &log_msg)
    {
        RHandleMoMALoggerMessage msg("[ERROR]", MoMALoggerLevel::ERRORS, MoMALogger::get_level());

        msg << log_msg;
        return msg;
    }

    static RHandleMoMALoggerMessage warning(const std::string &log_msg)
    {
        RHandleMoMALoggerMessage msg("[WARNING]", MoMALoggerLevel::WARNING,
                                     MoMALogger::get_level());

        msg << log_msg;
        return msg;
    }

    static RHandleMoMALoggerMessage message(const std::string &log_msg)
    {
        RHandleMoMALoggerMessage msg("[MESSAGE]", MoMALoggerLevel::MESSAGES,
                                     MoMALogger::get_level());
        msg << log_msg;
        return msg;
    }

    static MoMALoggerMessage info(const std::string &log_msg)
    {
        MoMALoggerMessage msg("[INFO]", MoMALoggerLevel::INFO, MoMALogger::get_level(),
                              MoMALogger::get_ostream());

        msg << log_msg;
        return msg;
    }

    static MoMALoggerMessage debug(const std::string &log_msg)
    {
        MoMALoggerMessage msg("[DEBUG]", MoMALoggerLevel::DEBUG, MoMALogger::get_level(),
                              MoMALogger::get_ostream());

        msg << log_msg;
        return msg;
    }

    static void set_level(MoMALoggerLevel logger_level)
    {
        MoMALogger::getInstance().logger_level = logger_level;
    }

    static MoMALoggerLevel get_level() { return MoMALogger::getInstance().logger_level; }

    static std::ostream &get_ostream() { return MoMALogger::getInstance().logger_ostream; }

  private:
    MoMALoggerLevel logger_level = MoMALoggerLevel::MESSAGES;
    std::ostream &logger_ostream = Rcpp::Rcout;

    // Default constructor ==> logger level = MESSAGES, output = Rcpp::Rcout

    // MESSAGES are things that the user should know, but doesn't need to
    // be concerned about => show them by default. Power users can suppress
    // them if desired.
    MoMALogger() {}

    MoMALogger(MoMALogger &&) = default;
    MoMALogger &operator=(MoMALogger &&) = default;
};

#endif
