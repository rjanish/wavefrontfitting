#include "abr_fit_image.h"
#include <cmath>
using namespace std;

void grad_frac_below_circle_wrt_dfposition(double x, double y, double R,
                                           double frac_below_grad[2], 
                                           double l)
{
    // initialize derivitive of fractional area to case of pixel above 
    // or below circle
    frac_below_grad[0] = 0.0;
    frac_below_grad[1] = 0.0;

    // 0 = x = y:
    if(y == 0)
    {
        frac_below_grad[0] = 0.0;
        frac_below_grad[1] = 0.0;
        return;
    }
        
    // for 0 <= x < y:  
    const double r = sqrt(pow(x, 2) + pow(y, 2));
    const double a = 0.5 + x/(2*y) + ((r - R)/l)*(r/y);
    const double DaDx = (0.5 + (x/l)*(2 - R/r))/y;
    const double DaDy = ((R - r)*r/l - x/2)/pow(y, 2) + (2 - R/r)/l;
        
    if ((a <= 0) or (a >= 1 + x/y))  // Case 1, pixel below circle, and 
                                     // Case 5, pixel above circle
    {
        frac_below_grad[0] = 0.0;
        frac_below_grad[1] = 0.0;
    }
    else if ((x/y <= a) and (a <= 1)) // Case 3, circle intersects left 
                                      // and right sides
    {
        frac_below_grad[0] = 1/(2*y) - DaDx;
        frac_below_grad[1] = -x/(2*pow(y, 2)) - DaDy;
    }
    else if ((0 < a) and (a < x/y)) // Case 2, circle intersects top 
                                    // and right
    {
        frac_below_grad[0] = -(y*a/x)*(DaDx - a/(2*x));
        frac_below_grad[1] = -(y*a/x)*(DaDy + a/(2*y));
    }
    else // Case 4: (1 < a) and (a < 1 + x/y), 
         // circle intersects left and bottom
    {
        frac_below_grad[0] = (1 - pow((a - 1)*y/x, 2))/(2*y) - \
                    (1 - (a - 1)*y/x)*DaDx;
        frac_below_grad[1] = -(1 - pow((a - 1)*y/x, 2))*x/(2*pow(y, 2)) - \
                    (1 - (a - 1)*y/x)*DaDy;
    }
    
    return;
}

double partial_frac_below_circle_wrt_defocus(double x, double y, double R, 
                                             double l, double defocus)
{
    // initialize derivitive of fractional area to case of pixel 
    // above or below circle
     double Df = 0;
    // for y = x = 0:
    if(y == 0)
    {
        Df = 0.0;
        return Df;
    }
    // for 0 <= x < y:
    double r = sqrt(pow(x,2) + pow(y,2));
    double a = 0.5 + x/(2*y) + ((r - R)/l)*(r/y);
    if ((a <= 0) or (a >= 1 + x/y)) // Case 1, pixel below circle, and 
                                    // Case 5, pixel above circle
    {
        Df = 0;
    }
    else if ((x/y <= a) and (a <= 1)) // Case 3, circle intersects left 
                                      // and right sides, x=0 calls this
    {
        Df = r*R/(l*y*defocus);
    }
    else if ((0 < a) and (a < x/y)) // Case 2, circle intersects top 
                                    // and right
    {
        Df = a*r*R/(l*x*defocus);
    }
    else // Case 4: (1 < a) and (a < 1 + x/y), 
         // circle intersects left and bottom
    {
        Df = (1 - (a - 1)*y/x)*r*R/(l*y*defocus);
    }
    
    return Df;
}


void grad_edge_term(double x_df, double y_df, double grad_x_df[],
                    double grad_y_df[], double params[], double l,
                    double Rin, double Rout, double edge_term[])
{
    // parameter lables
    const double defocus = params[7];
    // transformation flags
    int across_x = 0;
    int across_y = 0;
    int across_xequalsy = 0;
    // reflect x, y into pi/4 <= theta <= pi/2 region of pupil plane
    if (x_df < 0)
    {
        x_df *= -1;
        across_x = 1;
    }
    if (y_df < 0)
    {
        y_df *= -1;
        across_y = 1;
    }
    if (y_df < x_df)
    {
        const double tmp_holder = x_df;
        x_df = y_df;
        y_df = tmp_holder;
        across_xequalsy = 1;
    }

    // compute derivatives with repect to paramters 1 through 7
    for (int i = 1; i <= 7; i++)
    {

        // compute gradient of illuminated fraction using transformed 
        // coordinates with respect to x_dy, y_df and compute partial of 
        // illuminated fraction with respect to the defoucs parameter.
        double ill_frac_position_partials[2];
        double ill_frac_defocus_partial;
        const double r = sqrt(pow(x_df, 2) + pow(y_df, 2));
            // pixel near Rin
        if ((Rin - l/sqrt(2) < r) and (r < Rin + l/sqrt(2)))
        {
            grad_frac_below_circle_wrt_dfposition(x_df, y_df, Rin,
                                                  ill_frac_position_partials,
                                                  l);
            // include defoucs partial                             
            if (i == 7)
            {
                ill_frac_defocus_partial = \ 
                    partial_frac_below_circle_wrt_defocus(x_df, y_df, 
                                                          Rin, l, defocus);
            }
            else
            {
                ill_frac_defocus_partial = 0.0;                
            }
            // negate to get derivative of fraction outside circle
            for (int j = 0; j < 2; j++)
            {
                ill_frac_position_partials[j] *= -1;
            }
            ill_frac_defocus_partial *= -1;
        }
            // pixel near Rout
        else if ((Rout - l/sqrt(2) < r) and (r < Rout + l/sqrt(2)))
        {
            grad_frac_below_circle_wrt_dfposition(x_df, y_df, Rout,
                                                  ill_frac_position_partials,
                                                  l);
            // include defoucs partial                             
            if (i == 7)
            {
                ill_frac_defocus_partial = \ 
                    partial_frac_below_circle_wrt_defocus(x_df, y_df, 
                                                          Rout, l, defocus);
            }
            else
            {
                ill_frac_defocus_partial = 0.0;                
            }
        }
            // pixel does not intersect edge of illuminated area
        else
        {
            for (int j = 0; j < 2; j++)
            {
                ill_frac_position_partials[j] = 0.0;
            }
            ill_frac_defocus_partial = 0.0;
        }

        // Transform x and y derivatives back to origional coordinates
            // if y/x is not initially > 1, permute x and y
        if (across_xequalsy)
        {
            const double temp_var = ill_frac_position_partials[0];
            ill_frac_position_partials[0] = ill_frac_position_partials[1];
            ill_frac_position_partials[1] = temp_var;
        }
            // if coords not initially in quadrant I, reflect derivatives
            // back to initial quadrants
        if ((across_x == 1) and (across_y == 1))
        {
            ill_frac_position_partials[0] *= -1;
            ill_frac_position_partials[1] *= -1;
        }
        else if ((across_x == 1) and (across_y == 0))
        {
            ill_frac_position_partials[0] *= -1;
        }
        else if ((across_x == 0) and (across_y == 1))
        {
            ill_frac_position_partials[1] *= -1;
        }

        // total derivative of illuminated fraction with respect to 
        // parameter i via chain rule:
        edge_term[i] = ill_frac_position_partials[0]*grad_x_df[i] + \
                       ill_frac_position_partials[1]*grad_y_df[i] + \
                       ill_frac_defocus_partial;
    }

    // set derivative with respect to pupil_flux and background level to zero:
    edge_term[8] = 0.0;
    edge_term[9] = 0.0;

    return;
}
