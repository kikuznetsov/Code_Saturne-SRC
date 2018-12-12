#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_from_rije_profile.h"


#ifdef __cplusplus
extern "C" {
#endif
    
/**
* Returns the interpolated y-value.
* Saturates to y0 or y1 if x outside interval [x0, x1].
*/
struct record_rijssg_t
interpolate_rijssg_record( double y0,
                   struct record_rijssg_t* rec0,
                   double y1,
                  struct record_rijssg_t* rec1,
                  double y)
{
    double t;

    if (y <= y0) { return *rec0; }
    if (y >= y1) { return *rec1; }

    t =  (y-y0);
    t /= (y1-y0);

    struct record_rijssg_t temp;
    temp.y = y;
    temp.u = rec0->u + t*(rec1->u - rec0->u);
    temp.v = rec0->v + t*(rec1->v - rec0->v);

//    printf("rec0.u = %f\n",rec0->u);
//    printf("rec1.u = %f\n",rec1->u);
//    printf("temp.u = %f\n",temp.u);
    temp.rxx = rec0->rxx + t*(rec1->rxx-rec0->rxx);
    temp.ryy = rec0->ryy + t*(rec1->ryy-rec0->ryy);
    temp.rzz = rec0->rzz + t*(rec1->rzz-rec0->rzz);

    temp.rxy = rec0->rxy + t*(rec1->rxy-rec0->rxy);
    temp.ryz = rec0->ryz + t*(rec1->ryz-rec0->ryz);
    temp.rxz = rec0->rxz + t*(rec1->rxz-rec0->rxz);

    temp.eps = rec0->eps + t*(rec1->eps-rec0->eps);
    return temp;
}
/******************************************************************************/




struct record_rijssg_t
interpolate_rijssg(struct profile_rijssg_t* rows, double y_new)
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
            return interpolate_rijssg_record(rows->rec[segment].y,   /* x0 */
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

//int read_profile_SSG(const char *fName, size_t num_lines, struct profile_rijssg_t* rows) {
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
//        rows->rec[i].rxx = atof(rowFields[6]);
//        rows->rec[i].ryy = atof(rowFields[7]);
//        rows->rec[i].rzz = atof(rowFields[8]);
//        rows->rec[i].rxy = atof(rowFields[9]);
//        rows->rec[i].ryz = atof(rowFields[10]);
//        rows->rec[i].rxz = atof(rowFields[11]);
//        rows->rec[i].eps = atof(rowFields[12]);
//        i++;
//        //destroy line
//        CsvParser_destroy_row(row);
//        if (++line_count >= num_lines)
//          break;
//    }
//    CsvParser_destroy(csvparser);
//    return EXIT_SUCCESS;
//}


int read_profile_SSG(const char *fName, size_t num_lines, struct profile_rijssg_t* rows) {
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

        tok = strtok(NULL, ",");                //"rxx"
        rows->rec[i].rxx = atof(tok);

        tok = strtok(NULL, ",");                //"ryy"
        rows->rec[i].ryy = atof(tok);

        tok = strtok(NULL, ",");                //"rzz"
        rows->rec[i].rzz = atof(tok);

        tok = strtok(NULL, ",");                //"rxy"
        rows->rec[i].rxy = atof(tok);

        tok = strtok(NULL, ",");                //"ryz"
        rows->rec[i].ryz = atof(tok);

        tok = strtok(NULL, ",");                //"rxz"
        rows->rec[i].rxz = atof(tok);

        tok = strtok(NULL, ",");                //"eps"
        rows->rec[i].eps = atof(tok);

        i++;
        if (++line_count >= num_lines)
          break;
    }
    fclose(stream);
    return EXIT_SUCCESS;
}


#ifdef __cplusplus
}
#endif
