/*
http://suacommunity.com/dictionary/gettimeofday-entry.php
*/

#include <windows.h>
#include <time.h>

struct timezone
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz);