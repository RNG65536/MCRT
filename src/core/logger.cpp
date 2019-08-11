#include "logger.h"

static Logger s_debug_logger(false);
static Logger s_info_logger(true);

Logger::Logger(bool enabled) : m_enabled(enabled)
{
}

Logger& Logger::debug()
{
    return s_debug_logger;
}

Logger& Logger::info()
{
    return s_info_logger;
}
