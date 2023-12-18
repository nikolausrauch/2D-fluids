#pragma once

#include <chrono>
#include <map>
#include <string>

typedef std::conditional<std::chrono::high_resolution_clock::is_steady,
                         std::chrono::high_resolution_clock,
                         std::chrono::steady_clock >::type  HighResClock;
typedef HighResClock::time_point    TimePoint;
typedef HighResClock::duration      Time;

#define profile_sample(name) CPU_Sample sample_ ## name(#name)

class PerfMonitor
{
public:
    static PerfMonitor& instance();

    void start_frame();
    void clear();

    void add_time(const std::string& name, double time);
    double time(const std::string& name);

private:
    std::map<std::string, double> mTimes;
};

struct CPU_Sample
{
    CPU_Sample(const std::string& name);
    ~CPU_Sample();

private:
    const std::string mName;
    TimePoint mStart;
};
