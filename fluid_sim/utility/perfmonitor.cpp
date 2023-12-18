#include "perfmonitor.h"
#include <cmath>


PerfMonitor& PerfMonitor::instance()
{
    static PerfMonitor s_monitor;
    return s_monitor;
}

void PerfMonitor::start_frame()
{
    for(auto& times : mTimes)
    {
        times.second = 0;
    }
}

void PerfMonitor::clear()
{
    mTimes.clear();
}

void PerfMonitor::add_time(const std::string& name, double time)
{
    mTimes[name] += time;
}

double PerfMonitor::time(const std::string& name)
{
    return mTimes[name];
}


CPU_Sample::CPU_Sample(const std::string& name)
    : mName(name), mStart(HighResClock::now())
{

}

CPU_Sample::~CPU_Sample()
{
    auto end = HighResClock::now();
    PerfMonitor::instance().add_time(mName, std::chrono::duration<double, std::ratio<1, 1000>>(end - mStart).count());
}
