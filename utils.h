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

// this function reads a cvs file and modifies the array data, it works, do not touch it !
//the file sould has 3 rows
int read_csv(const char *filename, gsl_matrix* data);

// this function suppress the lines with zero detection in the csv dataset
gsl_matrix* suppress_zero_lines(const gsl_matrix* data);

// this function returns the rotation matrix following the axis 'axis' and the angle q
gsl_matrix* Rotd_axis(const int axis, const double q);

// this function returns the rotation matrix following the angles thetad,phid,psid
gsl_matrix* Rot_zxz(const double thetad, const double phid, const double psid);

// this function is the function we want to minimize it returns the standart deviation of the pixel offset
double cost_function(gsl_vector* X, void* parameters);

// this function will find the optimal camera parameters
int optimizer(Global_parameters parameters);


#endif //ROBOT_CAMERA_CALIBRATION_UTILS_H
