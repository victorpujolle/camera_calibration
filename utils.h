/*GSL*/
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>

/*Output to console*/
#include <iostream>

/*Computational time*/
#include <ctime>

/*The general headers*/
#include <pwd.h>
#include <stdio.h>
#include <iostream>
#include <signal.h>
#include <string.h>
#include <math.h>

/*Specific headers*/
#include "Global_parameters.h"

/*Special definitions*/
#define PI 3.14159265

using namespace std;

#ifndef ROBOT_CAMERA_CALIBRATION_UTILS_H
#define ROBOT_CAMERA_CALIBRATION_UTILS_H

// this function prints a gsl matrix
int print_matrix(const gsl_matrix *mat);

// this function prints a gsl vector
int print_vector(const gsl_vector* vect);

// computes the mean of the vector
double vector_mean(const gsl_vector* u);

// this function modifies the input vector u <- (u- mean(u)).^2
int vector_variance_element(gsl_vector* u);

// applies sqrt to each element of the vector u
int vector_sqrt_element(gsl_vector* u);

// this function reads a cvs file and modifies the array data, it works, do not touch it !
//the file sould has 3 rows
int read_csv(const char *filename, gsl_matrix* data);

// this function suppress the lines with zero detection in the csv dataset
gsl_matrix* suppress_zero_lines(const gsl_matrix* data);

// this function returns the rotation matrix following the axis 'axis' and the angle q
gsl_matrix* Rotd_axis(const int axis, const double q);

// this function returns the rotation matrix following the angles thetad,phid,psid
gsl_matrix* Rot_zxz(const double thetad, const double phid, const double psid);

// this function initialize the struct parameters
parameters set_parameters(gsl_matrix* calib_data, double range_dm=0.10, double range_dc=0.05, double range_df=30, double range_dpx=30, double range_angle=10, double inmarker_rotation=0, int minval=100, int repeat_no=50);

// this is the functions one should optimize
double _cost_function(gsl_vector* X, void* param);

// this function will find the optimal camera parameters
int optimizer(parameters* param);


#endif //ROBOT_CAMERA_CALIBRATION_UTILS_H
