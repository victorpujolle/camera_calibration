#include <pwd.h>
#include <stdio.h>
#include <iostream>

/*GSL*/
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

using namespace std;

#ifndef ROBOT_CAMERA_CALIBRATION_UTILS_STRUCT_H
#define ROBOT_CAMERA_CALIBRATION_UTILS_STRUCT_H

struct parameters
{
    gsl_matrix* calib_data;      // calibration data
    size_t len_data;
    double range_dm;        // deviation scale of mask position
    double range_dc;        // deviation scale of camera mount offset
    double range_df;        // deviation scale of focal length value
    double range_dpx;       // deviation scale of pixel offset
    double range_angle;     // deviation scale of pixel angular offsets
    double marker_rotation; // marker rotation - marker may positioned with 0, 90, 180, 270 degrees rotation
    int minval; // initial minimum value
    int repeat_no;
    gsl_matrix *px1,*px2, *M0, *f;    // values derived from transform matrix
    gsl_matrix* T_cb[]; // array of matrix

};

#endif //ROBOT_CAMERA_CALIBRATION_UTILS_STRUCT_H
