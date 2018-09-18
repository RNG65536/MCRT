#define NOMINMAX
#include <Windows.h>
#include "timer.h"

struct TimerContext
{
    LARGE_INTEGER hardware_freq, start_point, stop_point;
};

Timer::~Timer()
{
    delete m_ctx;
}

void Timer::reset(void)
{
    TimerContext* ctx = static_cast<TimerContext*>(m_ctx);
    QueryPerformanceCounter(&ctx->start_point);
    ctx->stop_point = ctx->start_point;
}

float Timer::getTime(void)
{
    TimerContext* ctx = static_cast<TimerContext*>(m_ctx);
    QueryPerformanceCounter(&ctx->stop_point);
    return float(ctx->stop_point.QuadPart - ctx->start_point.QuadPart) /
           float(ctx->hardware_freq.QuadPart);
}

Timer::Timer(void)
{
    auto ctx = new TimerContext();
    m_ctx = ctx;
    QueryPerformanceFrequency(&ctx->hardware_freq);
    QueryPerformanceCounter(&ctx->start_point);
    ctx->stop_point = ctx->start_point;
}
