#include <time.h>
#include <stdio.h>

class Timer {

 private:
  time_t start_time, stop_time;
  char message[80];

 public:

  Timer() { start_time = time(NULL); }
  ~Timer() { }
  
  // note: ctime will insert a \n at the end of it's string
  char* start() { start_time = time(NULL); sprintf(message, "Start Time: %s", ctime(&start_time)); return message; }
  char* stop() { stop_time = time(NULL); sprintf(message, "Stop Time: %s", ctime(&stop_time)); return message; }

  char* report() {
  
    time_t elapsed = (time_t) difftime(stop_time, start_time);
    struct tm* report = gmtime(&elapsed);
    sprintf(message, "Elapsed time: %d day(s), %d hour(s), %d min(s), %d sec(s)",
	    report->tm_yday, report->tm_hour, report->tm_min, report->tm_sec);

    return message;
  }
};

