// timer.h

#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

class Timer {

public:
  Timer();
  void Start();
  void Stop();
  double Duration();

  Clock::time_point start_time;
  Clock::time_point stop_time;
  bool is_stopped;
};

Timer::Timer() {
  Start();
}

void Timer::Start() {
  start_time = Clock::now();
  is_stopped = false;
}

void Timer::Stop() {
  stop_time = Clock::now();
  is_stopped = true;
}

double Timer::Duration() {
  std::chrono::duration<double> dur;
  if (is_stopped) {
    dur = stop_time - start_time;
  } else {
    Clock::time_point curr_time;
    curr_time = Clock::now();
    dur = curr_time - start_time;
  }
  return dur.count();
}
