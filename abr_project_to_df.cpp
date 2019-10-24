#include "abr_fit_image.h"
#include "nr_routines.h"
#include "nrutil.h"
#include <cmath>
using namespace std;

void grad_wavefront_residual(int n, double x[], double f[], double params[], 
                             double x_im, double y_im)
{
    // aberration parameters
    const double coma_x = params[1];
    const double coma_y = params[2];
    const double astig_x = params[3];
    const double astig_y = params[4];
    const double defocus = params[7];
    // x-component coefficents
    const double gradx_x_lin = 1.0 + astig_x/defocus;
    const double gradx_y_lin = astig_y/defocus;
    const double gradx_x_sq = 3*coma_x/(4*pow(defocus, 2));
    const double gradx_y_sq = coma_x/(4*pow(defocus, 2));
    const double gradx_cross = coma_y/(2*pow(defocus, 2));
    // y-component coefficents
    const double grady_x_lin = astig_y/defocus;
    const double grady_y_lin = 1.0 - astig_x/defocus;
    const double grady_x_sq = coma_y/(4*pow(defocus, 2));
    const double grady_y_sq = 3*coma_y/(4*pow(defocus, 2));
    const double grady_cross = coma_x/(2*pow(defocus, 2));
    // transfroms
    f[1] = (gradx_x_sq*x[1] + gradx_x_lin)*x[1] + \
           (gradx_y_sq*x[2] + gradx_y_lin + gradx_cross*x[1])*x[2] - x_im;

    f[2] = (grady_x_sq*x[1] + grady_x_lin)*x[1] + \
           (grady_y_sq*x[2] + grady_y_lin + grady_cross*x[1])*x[2] - y_im;
    return;
}


int project_to_defocus_only(double x_im, double y_im, double params[],
                            double *x_df, double *y_df)
{
    double *df_only_point = vector(1, 2);
    df_only_point[1] = x_im;
    df_only_point[2] = y_im;
    int convergence_flag = 0;
    double *f = vector(1, 2); // working space for newt routine
    const int dimension_of_eqn_system = 2; // magic number
    newt(df_only_point, dimension_of_eqn_system, &convergence_flag, params, 
         x_im, y_im, grad_wavefront_residual); // nr newton's method solver
    *x_df = df_only_point[1];
    *y_df = df_only_point[2]; 
    free_vector(f, 1, 2);
    free_vector(df_only_point, 1, 2);
    return convergence_flag;
}
