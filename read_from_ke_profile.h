#ifndef READ_FROM_KEPROFILE_H
#define READ_FROM_KEPROFILE_H

// #include <stdio.h>
// #include <stdlib.h>


#ifdef __cplusplus
extern "C" {
#endif


struct  record_keps_t{
    double y;
    double u;
    double v;
    double k;
    double eps;
};


struct profile_keps_t{
    size_t n_rows;
    struct record_keps_t* rec;
};


/**
* Returns the interpolated y-value.
* Saturates to y0 or y1 if x outside interval [x0, x1].
*/
struct record_keps_t
interpolate_keps_record( double y0,
                   struct record_keps_t* rec0,
                   double y1,
                  struct record_keps_t* rec1,
                  double y);

struct record_keps_t interpolate_keps(struct profile_keps_t* rows, double y_new);
int read_profile_keps(const char *fName, size_t num_lines, struct profile_keps_t* rows);

#ifdef __cplusplus
}
#endif

#endif // READ_FROM_KEPROFILE_H
