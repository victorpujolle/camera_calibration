#include "utils.h"


using namespace std;

// this function prints a gsl matrix
int print_matrix(const gsl_matrix *mat)
{
    if (mat == nullptr)
    {
        printf("matrix null\n");
        return 0;
    }

    size_t m = mat->size1;
    size_t n = mat->size2;

    for (int i = 0; i < m; i++)
    {
        printf("[ ");
        for (int j = 0; j < n; j++)
        {
            printf(" %.4f  ", gsl_matrix_get(mat,i,j));
        }
        printf(" ]\n");
    }
    printf("\n");
    return 0;
};

// this function prints a gsl vector
int print_vector(const gsl_vector* vect)
{
    if (vect == nullptr)
    {
        printf("vector null\n");
        return 0;
    }
    int m = vect->size;

    for (int i = 0; i < m; i++)
    {
        printf("[ %.4f ] ", gsl_vector_get(vect,i));
        printf("\n");
    }
    printf("\n");
    return 0;
};

// this function reads a cvs file and modifies the array data, it works, do not touch it !
//the file should has 3 rows
int read_csv(const char *filename, gsl_matrix* data)
{
    FILE* my_file = fopen(filename,"r");


    if (!my_file) { /* open operation failed. */
        perror("error ");
        printf("Failed opening file %s for reading",filename);
        exit(1);
    }

    char buffer[100];
    char* buffer2;

    int i = 0, j = 0;

    // this for read the csv file
    while (fscanf(my_file, "%s", buffer)!=EOF)
    {
        buffer2 = strtok (buffer,",");
        while (buffer2 != NULL)
        {
            gsl_matrix_set(data,i,j,atof(buffer2));
            buffer2 = strtok (NULL, ",");
            j++;
        }
        j = 0;
        i++;
    }
    fclose(my_file);

    return 0;


};

// this function suppress the lines with zero detection in the csv dataset
gsl_matrix* suppress_zero_lines(const gsl_matrix* data)
{
    size_t n = data->size1;
    size_t m = data->size2;
    double x;
    short int n_res = n;
    gsl_vector* zero_lines_index = gsl_vector_calloc(n);

    for(int i = 0 ; i < n ; i++)
    {
        for (int j = 0 ; j < 2 ; j++)
        {
            x = gsl_matrix_get(data,i,j);
            if(x != 0) break;
            if(j == 1)
            {
                gsl_vector_set(zero_lines_index,i,1);
                n_res--;
            }
        }
    }

    gsl_matrix* res = gsl_matrix_alloc(n_res,m);
    gsl_vector* line = gsl_vector_alloc(m);
    short int  index = 0;

    for (int i = 0 ; i < n ; i++)
    {
        if (gsl_vector_get(zero_lines_index, i) == 0)
        {
            gsl_matrix_get_row(line, data, i);
            gsl_matrix_set_row(res,index,line);

            index++;
        }
    }

    return res;
};

// this function returns the rotation matrix followi+ng the axis 'axis' and the angle q
gsl_matrix* Rotd_axis(const int axis, const double q)
{
    // conversion of q to radians
    double qrad = q*PI/180.0;
    // computation of cos(q) and sin(q)
    double cosq = cos(qrad);
    double sinq = sin(qrad);

    gsl_matrix* R = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(R);

    switch (axis)
    {
        case 1:
            gsl_matrix_set(R,1,1,cosq);
            gsl_matrix_set(R,1,2,-sinq);
            gsl_matrix_set(R,2,1,sinq);
            gsl_matrix_set(R,2,2,cosq);
            break;
        case 2:
            gsl_matrix_set(R,0,0,cosq);
            gsl_matrix_set(R,0,2,sinq);
            gsl_matrix_set(R,2,0,-sinq);
            gsl_matrix_set(R,2,2,cosq);
            break;
        case 3:
            gsl_matrix_set(R,0,0,cosq);
            gsl_matrix_set(R,0,1,-sinq);
            gsl_matrix_set(R,1,0,sinq);
            gsl_matrix_set(R,1,1,cosq);
            break;
    }

    return R;
};

// this function returns the rotation matrix following the angles thetad,phid,psid in the order zxz
gsl_matrix* Rot_zxz(const double thetad, const double phid, const double psid)
{
    gsl_matrix* Rtheta = Rotd_axis(3,thetad);
    gsl_matrix* Rphi = Rotd_axis(1,phid);
    gsl_matrix* Rpsi = Rotd_axis(3,psid);

    gsl_matrix* tmp1 = gsl_matrix_alloc(3,3);
    gsl_matrix* R = gsl_matrix_alloc(3,3);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Rtheta,Rphi,0.0,tmp1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,tmp1,Rpsi,0.0,R);

    return R;


};

// this function returns the rotation matrix following the angles zd, yd, xd in the order zyx
gsl_matrix* Rot_zyx(const double zd, const double yd, const double xd)
{
    gsl_matrix* Rz = Rotd_axis(3,zd);
    gsl_matrix* Ry = Rotd_axis(2,yd);
    gsl_matrix* Rx = Rotd_axis(1,xd);

    gsl_matrix* tmp1 = gsl_matrix_alloc(3,3);
    gsl_matrix* R = gsl_matrix_alloc(3,3);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Rz,Ry,0.0,tmp1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,tmp1,Rx,0.0,R);

    return R;


};

Global_parameters::Global_parameters(gsl_matrix* incalib_data, double inrange_dm, double inrange_dc, double inrange_df, double inrange_dpx, double inrange_angle, double inmarker_rotation, int inminval, int inrepeat_no)
{
    /* initialization of the global parameters */
    range_dm = inrange_dm;
    range_dc = inrange_dc;
    range_df = inrange_df;
    range_dpx = inrange_dpx;
    range_angle = inrange_angle;
    marker_rotation = inmarker_rotation;
    minval = inminval;
    repeat_no = inrepeat_no;

    len_data = incalib_data->size1;
    calib_data = incalib_data;
    /* setting of px */
    px1 = gsl_matrix_alloc(len_data,4);
    px2 = gsl_matrix_alloc(len_data,4);
    gsl_vector* col = gsl_vector_alloc(len_data);
    for (int i = 0; i < 8 ; ++i)
    {
        gsl_matrix_get_col(col,incalib_data,i);
        if(i%2) gsl_matrix_set_col(px2,i/2,col);
        else gsl_matrix_set_col(px1,i/2,col);
    }

    /* setting of M0 */
    M0 = gsl_matrix_alloc(len_data,3);
    for (int i = 8; i < 11; ++i)
    {
        gsl_matrix_get_col(col,calib_data,i);
        gsl_matrix_set_col(M0,i-8,col);
    }

    // setting of f
    f = gsl_matrix_alloc(len_data,2);
    for (int i = 11; i < 13; ++i)
    {
        gsl_matrix_get_col(col,calib_data,i);
        gsl_matrix_set_col(f,i-11,col);
    }

    /* setting of F0 */
    gsl_matrix* F0 = gsl_matrix_alloc(len_data,3);
    for (int i = 13; i < 16; ++i)
    {
        gsl_matrix_get_col(col,calib_data,i);
        gsl_matrix_set_col(F0,i-13,col);
    }

    /* setting of eu */
    gsl_matrix* eu = gsl_matrix_alloc(len_data,3);
    for (int i = 16; i < 19; ++i)
    {
        gsl_matrix_get_col(col,calib_data,i);
        gsl_matrix_set_col(eu,i-16,col);
    }

    /* setting of dist */
    gsl_matrix* dist = gsl_matrix_alloc(len_data,1);
    gsl_matrix_get_col(col,calib_data,19);
    gsl_matrix_set_col(dist,0,col);

    /* setting T_bf, T_ff, T_fc, T_bc and T_cb */
    gsl_matrix* T_bf[len_data];
    gsl_matrix* T_ff[len_data];
    gsl_matrix* T_fc[len_data];
    gsl_matrix* T_bc[len_data];

    /* setting of the temp matrix */
    gsl_matrix* Rota;
    gsl_matrix* T_bfT_ff = gsl_matrix_alloc(4,4);
    gsl_matrix* Rt = gsl_matrix_alloc(3,3);
    gsl_matrix* P = gsl_matrix_alloc(3,1);
    gsl_matrix* RtP = gsl_matrix_alloc(3,1);

    for(int i_d = 0 ; i_d < len_data ; i_d++)
    {
        T_bf[i_d] = gsl_matrix_alloc(4,4);
        T_ff[i_d] = gsl_matrix_alloc(4,4);
        T_fc[i_d] = gsl_matrix_alloc(4,4);
        T_bc[i_d] = gsl_matrix_alloc(4,4);
        T_cb[i_d] = gsl_matrix_alloc(4,4);
        gsl_matrix_set_identity(T_bf[i_d]);
        gsl_matrix_set_identity(T_ff[i_d]);
        gsl_matrix_set_identity(T_fc[i_d]);
        gsl_matrix_set_identity(T_cb[i_d]);

        /* T_bf */
        for(int j = 0 ; j < 3 ; j++)
        {
            gsl_matrix_set(T_bf[i_d],j,3,gsl_matrix_get(F0,i_d,j));
        }

        /* T_ff */
        Rota = Rot_zxz(gsl_matrix_get(eu,i_d,0),gsl_matrix_get(eu,i_d,1),gsl_matrix_get(eu,i_d,2));
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3 ; ++j)
            {
                gsl_matrix_set(T_ff[i_d],i,j,gsl_matrix_get(Rota,i,j));
            }
        }

        /* T_fc */
        gsl_matrix_set(T_fc[i_d],2,3,gsl_matrix_get(dist,i_d,0));

        /* T_bc */
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,T_bf[i_d],T_ff[i_d],0.0,T_bfT_ff);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,T_bfT_ff,T_fc[i_d],0.0,T_bc[i_d]);

        /* Rt */
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                gsl_matrix_set(Rt,i,j,gsl_matrix_get(T_bc[i_d],i,j));
            }
        }

        /* P */
        for (int i = 0; i < 3; ++i)
        {
            gsl_matrix_set(P,i,0,gsl_matrix_get(T_bc[i_d],i,3));
        }

        /* T_cb */
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                gsl_matrix_set(T_cb[i_d],i,j,gsl_matrix_get(Rt,j,i));
            }
        }

        gsl_blas_dgemm(CblasTrans,CblasNoTrans,-1.0,Rt,P,0.0,RtP);
        for (int i = 0; i < 3; ++i)
        {
            gsl_matrix_set(T_cb[i_d],i,3,gsl_matrix_get(RtP,i,0));
        }
    }
};

// this function is the function we want to minimize it returns the standart deviation of the pixel offset
double cost_function(gsl_vector* X, void* parameters)
{
    Global_parameters param = *(Global_parameters *) parameters;

    size_t len_data = param.len_data;

    /* multiply normalized calibration variables with scales */
    double dfx = gsl_vector_get(X,0)*param.range_df;     // focal length
    double dfy = gsl_vector_get(X,1)*param.range_df;
    double dcx = gsl_vector_get(X,2)*param.range_dpx;    // pixel offset
    double dcy = gsl_vector_get(X,3)*param.range_dpx;
    double dxc = gsl_vector_get(X,4)*param.range_dc;     // camera offset
    double dyc = gsl_vector_get(X,5)*param.range_dc;
    double dzc = gsl_vector_get(X,6)*param.range_dc;
    double dac = gsl_vector_get(X,7)*param.range_angle;  // camera angle offset
    double dbc = gsl_vector_get(X,8)*param.range_angle;
    double dcc = gsl_vector_get(X,9)*param.range_angle;
    double dxm = gsl_vector_get(X,10)*param.range_dm;    // marker offset
    double dym = gsl_vector_get(X,11)*param.range_dm;
    double dzm = gsl_vector_get(X,12)*param.range_dm;
    double dam = gsl_vector_get(X,13)*param.range_angle; // marker angle offset
    double dbm = gsl_vector_get(X,14)*param.range_angle;
    double dcm = gsl_vector_get(X,15)*param.range_angle;

    /* marker location */
    gsl_matrix* R_m = Rot_zxz(0,0,param.marker_rotation);

    /* camera offset */
    gsl_matrix* T_cc = gsl_matrix_calloc(4,4); // T_cc = [Rot_zxz(dac,dbc,dcc), [dxc;dyc;dzc];0,0,0,1]
    gsl_matrix* Rot_zxz_dac_dbc_dcc = Rot_zxz(dac,dbc,dcc);
    for (int i = 0 ; i < 3 ; i++)
    {
        for (int j = 0 ; j < 3 ; j++)
        {
            gsl_matrix_set(T_cc,i,j,gsl_matrix_get(Rot_zxz_dac_dbc_dcc,i,j));
        }
    }
    gsl_matrix_set(T_cc,0,3,dxc);
    gsl_matrix_set(T_cc,1,3,dyc);
    gsl_matrix_set(T_cc,2,3,dzc);
    gsl_matrix_set(T_cc,3,3,1);

    /* base offset */
    gsl_matrix* T_bb = gsl_matrix_calloc(4,4); // T_cc = [Rot_zxz(dac,dbc,dcc), [dxc;dyc;dzc];0,0,0,1]
    gsl_matrix* Rot_zyx_dam_dbm_dcm = Rot_zyx(dam,dbm,dcm);
    for (int i = 0 ; i < 3 ; i++)
    {
        for (int j = 0 ; j < 3 ; j++)
        {
            gsl_matrix_set(T_bb,i,j,gsl_matrix_get(Rot_zyx_dam_dbm_dcm,i,j));
        }
    }
    gsl_matrix_set(T_bb,0,3,dxm);
    gsl_matrix_set(T_bb,1,3,dym);
    gsl_matrix_set(T_bb,2,3,dzm);
    gsl_matrix_set(T_bb,3,3,1);

    /* corners: marker is 10cm square shape */
    gsl_matrix* corners = gsl_matrix_calloc(4,3);
    gsl_matrix_set(corners,0,0,-0.05); gsl_matrix_set(corners,0,1,0.05);
    gsl_matrix_set(corners,1,0,0.05);  gsl_matrix_set(corners,1,1,0.05);
    gsl_matrix_set(corners,2,0,0.05);  gsl_matrix_set(corners,2,1,-0.05);
    gsl_matrix_set(corners,3,0,-0.05); gsl_matrix_set(corners,3,1,-0.05);

    /* camera matrix */
    gsl_matrix* K_c = gsl_matrix_calloc(3,4);
    gsl_matrix_set(K_c,0,0, gsl_matrix_get(param.f,0,0)+dfx );
    gsl_matrix_set(K_c,0,2, 320+dcx );
    gsl_matrix_set(K_c,1,1, gsl_matrix_get(param.f,1,0)+dfy );
    gsl_matrix_set(K_c,1,2, 240+dcy );
    gsl_matrix_set(K_c,2,2, 1 );

    /* initializing the variables used during the iterations */
    gsl_vector* corners_j = gsl_vector_alloc(3); // corners(J,:)
    gsl_vector* M0_i = gsl_vector_alloc(3); // M0(I,:) first, then (M0(I,:)'+R_m*corners(J,:)')
    gsl_vector* P_cc = gsl_vector_alloc(4);
    gsl_vector* tmpP_cc = gsl_vector_alloc(4); // T_cc*squeeze(T_cb(I,:,:))*T_bb*[(M0(I,:)'+R_m*corners(J,:)');1]

    //for (int i = 0 ; i < len_data ; i++)
    for (int i = 0 ; i < 1 ; i++) // data index
    {
        for (int j = 0 ; j < 4 ; ++j) // marker corner position
        {
            gsl_matrix_get_row(corners_j,corners,j);
            gsl_matrix_get_row(M0_i,param.M0,i);
            gsl_blas_dgemv(CblasNoTrans,1.0,R_m,corners_j,1.0,M0_i);
            for (int k = 0 ; k < 3 ; k++)
            {
                gsl_vector_set(tmpP_cc,k,gsl_vector_get(M0_i,k));
            }
            gsl_vector_set(tmpP_cc,3,1.0);
            gsl_blas_dgemv(CblasNoTrans,1.0,T_bb,tmpP_cc,0.0,P_cc);
            gsl_blas_dgemv(CblasNoTrans,1.0,param.T_cb[i],P_cc,0.0,tmpP_cc);
            gsl_blas_dgemv(CblasNoTrans,1.0,T_cc,tmpP_cc,0.0,P_cc);

            if(gsl_)
        }

    }



}

// this function will find the optimal camera parameters
int optimizer(Global_parameters parameters)
{
    /* setup of the random generator */
    gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus); // the generator is a Tausworthe generator
    gsl_rng_set(rand_gen,0); // the seed used is the default seed

    /* setup of the optimizer */
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; // the Nelder-Mead Simplex algorithm is used
    gsl_multimin_fminimizer* optimizer = gsl_multimin_fminimizer_alloc(T, 16);

    /* set the initial step size to 1.0 */
    gsl_vector* ss = gsl_vector_alloc (16);
    gsl_vector_set_all (ss, 1.0);


    /* initialization of the starting point X0 */
    const short int sizeX0 = 16;
    gsl_vector* X0 = gsl_vector_alloc(sizeX0);
    double random_double;

    /* loop because the algorithm used is the Nelder-Mead method which is really sensitive to local minimum */
    //for (int i = 0 ; i < parameters.repeat_no ; i++)
    for (int i = 0 ; i < 1 ; i++)
    {
        /* X0 is initialized randomly for each iteration */
        for(int k = 0 ; k < sizeX0 ; k++)
        {
            gsl_vector_set(X0,k,gsl_rng_uniform(rand_gen)-0.5);
        }
        parameters.marker_rotation = 180;
        gsl_vector_set_all(X0,0.5);
        cost_function(X0,&parameters);



    }

    return 0;
}
