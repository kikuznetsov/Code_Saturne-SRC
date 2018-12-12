#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_from_ke_profile.h"


#ifdef __cplusplus
extern "C" {
#endif

/**
* Returns the interpolated y-value.
* Saturates to y0 or y1 if x outside interval [x0, x1].
*/
struct record_keps_t
interpolate_keps_record( double y0,
                   struct record_keps_t* rec0,
                   double y1,
                  struct record_keps_t* rec1,
                  double y)
{
    double t;

    if (y <= y0) { return *rec0; }
    if (y >= y1) { return *rec1; }

    t =  (y-y0);
    t /= (y1-y0);

    struct record_keps_t temp;
    temp.y = y;
    temp.u = rec0->u + t*(rec1->u - rec0->u);
    temp.v = rec0->v + t*(rec1->v - rec0->v);

//    printf("rec0.u = %f\n",rec0->u);
//    printf("rec1.u = %f\n",rec1->u);
//    printf("temp.u = %f\n",temp.u);
    temp.k = rec0->k + t*(rec1->k-rec0->k);
    temp.eps = rec0->eps + t*(rec1->eps-rec0->eps);
    return temp;
}
/******************************************************************************/




struct record_keps_t
interpolate_keps(struct profile_keps_t* rows, double y_new)
/* 1D Table lookup with interpolation */
{
    size_t segment;
    size_t len = rows->n_rows;
    /* Check input bounds and saturate if out-of-bounds */
    if (y_new > (rows->rec[len-1].y)) {
       /* y-value too large, saturate to max y-value */
        return rows->rec[len-1];
    }
    else if (y_new < (rows->rec[0].y)) {
       /* x-value too small, saturate to min y-value */
        return rows->rec[0];
    }

    /* Find the segment that holds x */
    for (segment = 0; segment<(len-1); segment++)
    {
        if ((rows->rec[segment].y   <= y_new) &&
            (rows->rec[segment+1].y >= y_new))
        {
            /* Found the correct segment */
            /* Interpolate */
//            printf("%f < y_new < %f\n",rows->rec[segment].y,rows->rec[segment+1].y);
//            printf("%f < u_new < %f\n",rows->rec[segment].u,rows->rec[segment+1].u);
////            struct record_rijssg_t* seg0 = &(rows->rec[segment]);
////            struct record_rijssg_t* seg1 = &(rows->rec[segment + 1]);
            return interpolate_keps_record(rows->rec[segment].y,   /* x0 */
                                      &rows->rec[segment],   /* y0 */
                                      rows->rec[segment+1].y, /* x1 */
                                      &rows->rec[segment+1], /* y1 */
                                      y_new);                         /* x  */
        }
    }

    /* Something with the data was wrong if we get here */
    /* Saturate to the max value */
    return rows->rec[len-1];
}

int read_profile_keps(const char *fName, size_t num_lines, struct profile_keps_t* rows) {
    FILE* stream = fopen(fName, "r");
    if (stream == NULL)
            return EXIT_FAILURE;
    char line[1024];

    const char *tok;
    size_t line_count = 0;
    size_t i=0;
    fgets(line, 1024, stream);          //read header
    while (fgets(line, 1024, stream) != NULL){
        //put everything in profile
        tok = strtok(line, ",");                //"s"
        //skip

        tok = strtok(NULL, ",");                //"x"
        //skip

        tok = strtok(NULL, ",");                //"y"
        rows->rec[i].y = atof(tok);

        tok = strtok(NULL, ",");                //"z"
        //skip

        tok = strtok(NULL, ",");                //"u"
        rows->rec[i].u = atof(tok);

        tok = strtok(NULL, ",");                //"v"
        rows->rec[i].v = atof(tok);

        tok = strtok(NULL, ",");                //"k"
        rows->rec[i].k = atof(tok);

        tok = strtok(NULL, ",");                //"eps"
        rows->rec[i].eps = atof(tok);

        i++;
        if (++line_count >= num_lines)
          break;
    }
    fclose(stream);
    return EXIT_SUCCESS;

}

//int read_profile_keps(const char *fName, size_t num_lines, struct profile_keps_t* rows) {
//    CsvParser *csvparser = CsvParser_new(fName, ",", 1);
//    CsvRow *header;
//    CsvRow *row;
//    size_t line_count = 0;
//    size_t i=0;

//    while ((row = CsvParser_getRow(csvparser)) ) {
//        char **rowFields = CsvParser_getFields(row);
//        //put everything in profile
//        rows->rec[i].y = atof(rowFields[2]);
//        rows->rec[i].u = atof(rowFields[4]);
//        rows->rec[i].v = atof(rowFields[5]);
//        rows->rec[i].k = atof(rowFields[6]);
//        rows->rec[i].eps = atof(rowFields[7]);
//        i++;
//        //destroy line
//        CsvParser_destroy_row(row);
//        if (++line_count >= num_lines)
//          break;
//    }
//    CsvParser_destroy(csvparser);
//    return EXIT_SUCCESS;
//}

#ifdef __cplusplus
}
#endif
