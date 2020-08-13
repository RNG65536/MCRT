#include <mutex>
#include "logger.h"

Logger::Logger(bool enabled, nanolog::LogLevel log_level)
    : m_enabled(enabled), m_log_level(log_level)
{
    if (m_enabled)
    {
        static std::once_flag flag;
        std::call_once(flag, []() {
            nanolog::initialize(
                nanolog::GuaranteedLogger(), "./", "log.txt", 1);
        });

        m_logger = std::make_shared<nanolog::NanoLogLine>(
            m_log_level, __FILE__, __func__, __LINE__);
    }
}

Logger Logger::debug()
{
    return Logger(false, nanolog::LogLevel::INFO);
}

Logger Logger::info()
{
    return Logger(true, nanolog::LogLevel::INFO);
}
