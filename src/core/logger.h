#pragma once

#include <NanoLog.hpp>
#include <iostream>
#include <memory>
#include <sstream>
#include <functional>

// for console output, use std::endl
// for logger output, endl is automatically handled
class Logger
{
public:
    static std::ostream& carriage_return(std::ostream& os)
    {
        os << '\r';
        return os;
    }

    Logger(bool              enabled = true,
           nanolog::LogLevel log_level = nanolog::LogLevel::INFO);

    template <typename T>
    Logger& operator<<(const T& t)
    {
        if (m_enabled)
        {
            std::cout << t;
            if (nanolog::is_logged(m_log_level))
            {
                (*m_logger) << t;
            }
        }
        return *this;
    }

    Logger& operator<<(std::ostream& (*f)(std::ostream&))
    {
        if (m_enabled)
        {
            f(std::cout);
            //if (nanolog::is_logged(m_log_level))
            //{
            //    std::stringstream ss;
            //    f(ss);
            //    (*m_logger) << ss.str();
            //}
        }
        return *this;
    }

    ~Logger()
    {
        if (m_enabled)
        {
            if (nanolog::is_logged(m_log_level))
            {
                // finish the line
                nanolog::NanoLog() == (*m_logger);
            }
        }
    }

    static Logger debug();
    static Logger info();

private:
    bool                                  m_enabled = true;
    std::shared_ptr<nanolog::NanoLogLine> m_logger;
    nanolog::LogLevel                     m_log_level;
};
