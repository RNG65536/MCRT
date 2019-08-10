#pragma once

#include <iostream>

namespace debug
{
class Logger
{
public:
    Logger(bool enabled = true);

    // this does not work
    template <typename T>
    Logger& operator<<(const T& t)
    {
        if (m_enabled)
        {
            std::cout << t;
        }
        return *this;
    }

    Logger& operator<<(std::ostream& (*f)(std::ostream&))
    {
        if (m_enabled)
        {
            f(std::cout);
        }
        return *this;
    }

    static Logger& debug();
    static Logger& info();

private:
    bool m_enabled = true;
};
}
