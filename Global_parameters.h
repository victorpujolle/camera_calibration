#include <pwd.h>
#include <stdio.h>
#include <iostream>

/*GSL*/
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

using namespace std;

#ifndef ROBOT_CAMERA_CALIBRATION_UTILS_STRUCT_H
#define ROBOT_CAMERA_CALIBRATION_UTILS_STRUCT_H

// this struct iswhere the global parameters are stored
class Global_parameters
{
public:
    gsl_matrix* calib_data;      // calibration data
    size_t len_data;
    double range_dm;        // deviation scale of mask position
    double range_dc;        // deviation scale of camera mount offset
    double range_df;        // deviation scale of focal length value
    double range_dpx;       // deviation scale of pixel offset
    double range_angle;     // deviation scale of pixel angular offsets
    double marker_rotation; // marker rotation - marker may positioned with 0, 90, 180, 270 degrees rotation
    gsl_matrix *px1,*px2, *M0, *f;    // values derived from transform matrix
    gsl_matrix* T_cb[]; // array of matrix
    int minval; // initial minimum value
    int repeat_no;
                                                         
    Global_parameters(gsl_matrix* calib_data, double range_dm=0.10, double range_dc=0.05, double range_df=30, double range_dpx=30, double range_angle=10, double inmarker_rotation=0, int minval=100, int repeat_no=50);
    double cost_function(gsl_vector* X);
    int optimizer();
};

#endif //ROBOT_CAMERA_CALIBRATION_UTILS_STRUCT_H
