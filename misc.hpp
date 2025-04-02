#ifndef TPG_MISC_HPP
#define TPG_MISC_HPP

#include <sys/time.h>

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? a : b)
#endif

double realtime();

#endif
