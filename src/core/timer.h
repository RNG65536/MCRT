#pragma once

class Timer
{
public:
    Timer(void);

    // get elapsed time from last reset() or class construction. return The
    // elapsed time
    float getTime(void);

    // reset the timer
    void reset(void);

    ~Timer();

private:
    void* m_ctx;
};
