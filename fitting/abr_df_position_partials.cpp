#include "abr_fit_image.h"
#include <cmath>
using namespace std;

void grad_df_position(double x_df, double y_df, double params[], 
                      double Hinv[][2], double grad_x_df[], 
                      double grad_y_df[])
{
    // convinent labels
    const double coma_x = params[1];
    const double coma_y = params[2];
    const double astig_x = params[3];
    const double astig_y = params[4];
    const double defocus = params[7];

    // compute derivative of defocus only position with respect to
    // aberration parameters 1 thorugh 7. The derivative is given by the
    // inverse hessian matix of the wavefront times a vector whoes
    // functional form is dependent on which parameter the differentiation
    // is with respect to.
    for (int i = 1; i <= 7; i++)
    {
        // compute aberration paramter dependent vector
        double abr_vector[2];  
        if (i == 1) // wrt coma x
        {
            abr_vector[0] = (-1.0/(4*pow(defocus,2)))*(3.0*pow(x_df, 2) + \
                                                       pow(y_df, 2));
            abr_vector[1] = -x_df*y_df/(2*pow(defocus,2));
        }
        else if (i == 2) // wrt coma y
        {
            abr_vector[0] = -x_df*y_df/(2*pow(defocus,2));
            abr_vector[1] = (-1.0/(4*pow(defocus,2)))*(pow(x_df, 2) + \
                                                       3.0*pow(y_df, 2));
        }
        else if (i == 3) // wrt astig x
        {
            abr_vector[0] = -x_df/defocus;
            abr_vector[1] = y_df/defocus;
        }
        else if (i == 4) // wrt astig y
        {
            abr_vector[0] = -y_df/defocus;
            abr_vector[1] = -x_df/defocus;
        }
        else if (i == 5) // wrt tilt x
        {
            abr_vector[0] = -1;
            abr_vector[1] = 0;
        }
        else if (i == 6) // wrt tilt y
        {
            abr_vector[0] = 0;
            abr_vector[1] = -1;
        }
        else if (i == 7) // wrt defocus
        {
            abr_vector[0] = 3*coma_x*pow(x_df,2)/(2*pow(defocus,3)) + \
                            coma_x*pow(y_df,2)/(2*pow(defocus,3)) + \
                            coma_y*x_df*y_df/pow(defocus,3) + \
                            astig_x*x_df/pow(defocus,2) + \
                            astig_y*y_df/pow(defocus,2); 
            abr_vector[1] = coma_y*pow(x_df,2)/(2*pow(defocus,3)) + \
                            3*coma_y*pow(y_df,2)/(2*pow(defocus,3)) + \
                            coma_x*x_df*y_df/pow(defocus,3) + \
                            astig_y*x_df/pow(defocus,2) - \
                            astig_x*y_df/pow(defocus,2);
        }
        // compute partials with respect to parameter i
        grad_x_df[i] = Hinv[0][0]*abr_vector[0] + Hinv[0][1]*abr_vector[1];
        grad_y_df[i] = Hinv[1][0]*abr_vector[0] + Hinv[1][1]*abr_vector[1];
    }

    // set derivatives wrt pupil_flux and background
    grad_x_df[8] = 0.0;
    grad_y_df[8] = 0.0;
    grad_x_df[9] = 0.0;
    grad_y_df[9] = 0.0;

    return;
}
