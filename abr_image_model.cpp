#include "abr_fit_image.h"
#include <cmath>
using namespace std;

#define R_IN 0.8251
#define R_OUT 1.85
#define NUM_SINGLE_UNBLUR_PARAMS 9

void unblurred_image_model(double x_im, double y_im, double abr_params[],
                           int num_abr_params, double *brightness, 
                           double partials[], double pixel_width)
{   

    // assign parameter names
    const double tilt_x = abr_params[5];
    const double tilt_y = abr_params[6];
    const double defocus = abr_params[7];
    const double pupil_flux = abr_params[8];
    const double background = abr_params[9];
        
    // translate image to account for tilts
    x_im = x_im - tilt_x;
    y_im = y_im - tilt_y;

    // Project to defocus-only, if no corresponding defocus-only point
    // exits, set pixel to background, zero everything else, and exit.
    int no_convergence = 0;
    double x_df = 0.0, y_df = 0.0;
    no_convergence = project_to_defocus_only(x_im, y_im, abr_params,
                                             &x_df, &y_df);
    if (no_convergence)
    {
        *brightness = background;
        for (int i = 1; i < num_abr_params; i++) // set partials wrt all params except the
        {                                        // background level to zero
            partials[i] = 0.0;
        }
        partials[9] = 1.0; // partial wrt background level
        return;
    }   

    // set radii of illuminated area in defocus-only image
    const double r_in_df = R_IN*2.0*abs(defocus);
    const double r_out_df = R_OUT*2.0*abs(defocus);

    // check for case of non-illuminated image pixel
    const double illuminated_fraction = ill_frac(x_df, y_df, r_in_df,
                                                 r_out_df, pixel_width);
    if (illuminated_fraction == 0)
    {
        *brightness = background;
        for (int i = 1; i < num_abr_params; i++) // set partial wrt all params except the 
        {                                        // background level to zero
            partials[i] = 0.0;
        }
        partials[9] = 1.0;  // partial wrt background level 
        return;
    }                                              

    // compute brightness
    const double Hinv_det = inverse_hessian_det(x_df, y_df, abr_params);
    const double df_only_flux = pupil_flux/(4.0*pow(defocus, 2));
    const double star_brightness = abs(Hinv_det)*illuminated_fraction*df_only_flux;
    *brightness = star_brightness + background;

    // compute derivatives
        // intermediate computations
    double H_inv[2][2];
    inverse_hessian(x_df, y_df, abr_params, H_inv, Hinv_det);
    double grad_x_df[NUM_SINGLE_UNBLUR_PARAMS + 1];
    double grad_y_df[NUM_SINGLE_UNBLUR_PARAMS + 1];
    grad_df_position(x_df, y_df, abr_params, H_inv, grad_x_df, grad_y_df);

        // bulk contrubution
    grad_bulk_term(x_df, y_df, grad_x_df, grad_y_df, 
                   H_inv, abr_params, partials);

        // if pixel is near edge, add in edge contrubution
    if (illuminated_fraction < 1)
    {
        double edge_term[NUM_SINGLE_UNBLUR_PARAMS + 1];
        grad_edge_term(x_df, y_df, grad_x_df, grad_y_df, abr_params, 
                       pixel_width, r_in_df, r_out_df, edge_term);  
        for (int i=1; i <= 7; i++) // pupil_flux and background components (8 and 9) 
        {                          // are not affected
            partials[i] += -edge_term[i]/illuminated_fraction;
        }
    }    

        // add in contrbution for defocus, pupil_flux and background parameters
    partials[7] += 2.0/defocus;
    partials[8] = -1.0/pupil_flux;
    partials[9] = 1.0;

        // scale all params except background by negative of star brightness
    for (int i=1; i < num_abr_params; i++)
    {
        partials[i] *= -star_brightness;
    }


    return;
}

#undef R_IN
#undef R_OUT
#undef NUM_SINGLE_UNBLUR_PARAMS
