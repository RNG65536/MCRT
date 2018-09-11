#include "timer.h"

void WTimer::init(void)
{
    QueryPerformanceFrequency(&_freq);
    QueryPerformanceCounter(&_start);
    _stop = _start;
}

WTimer::~WTimer()
{

}

void WTimer::reset(void)
{
    QueryPerformanceCounter(&_start);
    _stop = _start;
}

float WTimer::getTime(void)
{
    QueryPerformanceCounter(&_stop);
    return float(_stop.QuadPart - _start.QuadPart) / float(_freq.QuadPart);
}

WTimer::WTimer(void)
{
    init();
}
