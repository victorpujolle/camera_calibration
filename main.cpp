#include <gsl/gsl_matrix.h>
#include "utils.h"
#include "Global_parameters.h"

int main() {
    printf("- starting -- camera calibration --\n");

    char name_file[] = "../Calib_list/1.csv";
    gsl_matrix* file_data = gsl_matrix_alloc(22,20);
    read_csv(name_file , file_data);
    printf("- reading the file -- %s -- done\n",name_file);

    gsl_matrix* data = suppress_zero_lines(file_data);
    printf("- deleting zero detection lines done -- %lu lines deleted --\n",file_data->size1-data->size1);

    // setup the first global parameters of the problem
    //Global_parameters param = Global_parameters(data);
    //printf("- setting the global parameters done\n");

    // setup the first global parameters of the problem
    parameters param = set_parameters(data);
    printf("- setting the global parameters done\n");

    //Optimization of the parameters
    optimizer(&param);
    printf("- optimizing the parameters done\n");

    printf("- end\n");

    return 0;
}