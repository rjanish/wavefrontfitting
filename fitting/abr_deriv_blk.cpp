#include "abr_fit_image.h"
#include <cmath>
using namespace std;

void grad_bulk_term(double x_df, double y_df, double grad_x_df[],
                    double grad_y_df[], double Hinv[][2],
                    double params[], double bulk_term[])
{
    // convinent labels
    const double coma_x = params[1];
    const double coma_y = params[2];
    const double astig_x = params[3];
    const double astig_y = params[4];
    const double defocus = params[7];

    // Compute derivative of hessian with respect to all parameters 
    // except pupil_flux (param[1] to param[7])
    for (int i = 1; i <= 7; i++)
    {
        // Compute derivative of Hessian; this has two terms, one is
        // dependent on the parameter with respect to which we are
        // differentiating and the other is independing of this parameter
            // aberration parameter dependent part
        double DH_Dparam[2][2];
        if (i == 1) // coma_x
        {
            DH_Dparam[0][0] = 3.0*x_df/(2*pow(defocus,2));
            DH_Dparam[0][1] = y_df/(2*pow(defocus,2));
            DH_Dparam[1][0] = y_df/(2*pow(defocus,2));
            DH_Dparam[1][1] = x_df/(2*pow(defocus,2));
        }
        else if (i == 2) // coma_y
        {
            DH_Dparam[0][0] = y_df/(2*pow(defocus,2));
            DH_Dparam[0][1] = x_df/(2*pow(defocus,2));
            DH_Dparam[1][0] = x_df/(2*pow(defocus,2));
            DH_Dparam[1][1] = 3.0*y_df/(2*pow(defocus,2));
        }
        else if (i == 3) // astig_x
        {
            DH_Dparam[0][0] = 1.0/defocus;
            DH_Dparam[0][1] = 0.0;
            DH_Dparam[1][0] = 0.0;
            DH_Dparam[1][1] = -1.0/defocus;
        }
        else if (i == 4) // astig_y
        {
            DH_Dparam[0][0] = 0.0;
            DH_Dparam[0][1] = 1.0/defocus;
            DH_Dparam[1][0] = 1.0/defocus;
            DH_Dparam[1][1] = 0.0;
        }
        else if ((i == 5) or (i == 6)) // tilt_x or tilt_y
        {                                               
            DH_Dparam[0][0] = 0.0;
            DH_Dparam[0][1] = 0.0;
            DH_Dparam[1][0] = 0.0;
            DH_Dparam[1][1] = 0.0;
        }
        else if (i == 7) // defocus
        {
            DH_Dparam[0][0] = -3*coma_x*x_df/pow(defocus, 3) - \
                                coma_y*y_df/pow(defocus, 3) - \
                                astig_x/pow(defocus, 2);
            DH_Dparam[0][1] = -coma_y*x_df/pow(defocus, 3) - \
                                coma_x*y_df/pow(defocus, 3) - \
                                astig_y/pow(defocus, 2);
            DH_Dparam[1][0] = -coma_y*x_df/pow(defocus, 3) - \
                                coma_x*y_df/pow(defocus, 3) - \
                                astig_y/pow(defocus, 2);
            DH_Dparam[1][1] = -coma_x*x_df/pow(defocus, 3) - \
                                3*coma_y*y_df/pow(defocus, 3) + \
                                astig_x/pow(defocus, 2);
        }
            // add in aberration parameter independent part
        DH_Dparam[0][0] += (3.0*coma_x*grad_x_df[i] + \
                            coma_y*grad_y_df[i])/(2*pow(defocus,2)); 
        DH_Dparam[0][1] += (coma_y*grad_x_df[i] + \
                            coma_x*grad_y_df[i])/(2*pow(defocus,2));
        DH_Dparam[1][0] += (coma_y*grad_x_df[i] + \
                            coma_x*grad_y_df[i])/(2*pow(defocus,2));
        DH_Dparam[1][1] += (coma_x*grad_x_df[i] + \
                            3.0*coma_y*grad_y_df[i])/(2*pow(defocus,2));
        // compute: bulk term = trace(hessian^-1 * d_hessian/d_parameter)
        bulk_term[i] = Hinv[0][0]*DH_Dparam[0][0] + \
                       Hinv[0][1]*DH_Dparam[1][0] + \
                       Hinv[1][0]*DH_Dparam[0][1] + \
                       Hinv[1][1]*DH_Dparam[1][1];
    }
    // compute derivative with respect to flux and background paramter
    bulk_term[8] = 0.0;
    bulk_term[9] = 0.0;
    return;
}
