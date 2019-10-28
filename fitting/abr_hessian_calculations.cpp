#include "abr_fit_image.h"
#include <cmath>
using namespace std;

double inverse_hessian_det(double x_df, double y_df, double params[])
{

    // aberration parameters
    const double coma_x = params[1];
    const double coma_y = params[2];
    const double astig_x = params[3];
    const double astig_y = params[4];
    const double defocus = params[7];

    // compute global paramters for hessian determinant
    const double hd_constant = \
      1 - (pow(astig_x, 2) + pow(astig_y, 2))/pow(defocus, 2);
    const double hd_x_lin = \
      (1.0/pow(defocus, 2))*(coma_x*(2 - astig_x/defocus) - \
                             coma_y*astig_y/defocus);
    const double hd_y_lin = \
        (1.0/pow(defocus, 2))*(-coma_x*astig_y/defocus + \
                                coma_y*(2 + astig_x/defocus));
    const double hd_x_sq = \
        (1.0/(4*pow(defocus, 4)))*(3*pow(coma_x, 2) - pow(coma_y, 2));
    const double hd_y_sq = \
        (1.0/(4*pow(defocus, 4)))*(-pow(coma_x, 2) + 3*pow(coma_y, 2));
    const double hd_cross = (2.0/pow(defocus, 4))*coma_x*coma_y;
    // compute determinant of hessian matrix
    const double hd = (hd_x_sq*x_df + hd_x_lin)*x_df + \
                      (hd_y_sq*y_df + hd_y_lin + hd_cross*x_df)*y_df + \  
                       hd_constant; 
    // det(H^-1) = 1/det(H)
    return (1.0/(hd));
}

void inverse_hessian(double x_df, double y_df, double params[], 
                     double Hinv[][2], double Hinv_determinant)
{
    // aberration parameters
    const double coma_x = params[1];
    const double coma_y = params[2];
    const double astig_x = params[3];
    const double astig_y = params[4];
    const double defocus = params[7];

    Hinv[0][0] = Hinv_determinant*(coma_x*x_df/(2*pow(defocus,2)) + \
                                   3*coma_y*y_df/(2*pow(defocus,2)) + \
                                   1 - astig_x/defocus);
    Hinv[1][1] = Hinv_determinant*(3*coma_x*x_df/(2*pow(defocus,2)) + \
                                   coma_y*y_df/(2*pow(defocus,2)) + \
                                   1 + astig_x/defocus);
    Hinv[0][1] = (-1)*Hinv_determinant*(coma_x*y_df/(2*pow(defocus,2)) + \
                                        coma_y*x_df/(2*pow(defocus,2)) + \
                                        astig_y/defocus);
    Hinv[1][0] = Hinv[0][1];

    return;
}