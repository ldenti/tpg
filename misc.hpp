#ifndef PS_MISC_HPP
#define PS_MISC_HPP

#include <sys/time.h>

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? a : b)
#endif

double realtime();

#endif
