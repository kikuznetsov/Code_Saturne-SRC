#ifndef READ_FROM_PROFILE_H
#define READ_FROM_PROFILE_H


#ifdef __cplusplus
extern "C" {
#endif


struct record_rijssg_t{
    double y;
    double u;
    double v;
    double rxx;
    double ryy;
    double rzz;

    double rxy;
    double ryz;
    double rxz;

    double eps;
};

struct profile_rijssg_t{
    size_t n_rows;
    struct record_rijssg_t* rec;
};


/**
* Returns the interpolated y-value.
* Saturates to y0 or y1 if x outside interval [x0, x1].
*/
struct record_rijssg_t
interpolate_rijssg_record( double y0,
                   struct record_rijssg_t* rec0,
                   double y1,
                  struct record_rijssg_t* rec1,
                  double y);

struct record_rijssg_t
    interpolate_rijssg(struct profile_rijssg_t* rows, double y_new);

int read_profile_SSG(const char *fName, size_t num_lines, struct profile_rijssg_t* rows);
// int read_profile_keps(const char *fName, struct record_keps_t* rows);

#ifdef __cplusplus
}
#endif

#endif // READ_FROM_PROFILE_H
