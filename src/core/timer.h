#ifndef timer_h__
#define timer_h__

#define NOMINMAX
#include <Windows.h>

class WTimer
{
    LARGE_INTEGER _freq, _start, _stop;

public:
    WTimer(void);

    /* Get elapsed time from last reset()
    or class construction.
    return The elapsed time
    */

    float getTime(void);

    /*
    Reset the timer
    */

    void reset(void);

    ~WTimer();

private:
    void init(void);
};

#endif // timer_h__
