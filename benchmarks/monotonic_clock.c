/* Original benchmark helper: monotonic elapsed wall time for full R API calls. */
#define _POSIX_C_SOURCE 200809L
#include <time.h>
#include <R.h>
#include <Rinternals.h>

SEXP coloc_monotonic(void) {
    struct timespec value;
    if (clock_gettime(CLOCK_MONOTONIC, &value) != 0)
        error("clock_gettime(CLOCK_MONOTONIC) failed");
    return ScalarReal((double)value.tv_sec + (double)value.tv_nsec * 1e-9);
}
