#include "abr_fit_image.h"
#include <cmath>
using namespace std;

double frac_below_circle(double x, double y, double R, double l)
{   
    // This function is only valid if the x and y coodinates lie in the polar angle range pi/4 < theta < pi/2 and if R >> l

double frac_area = 0;
        // for y = x = 0:
    if(y == 0)
    {
        frac_area = 1.0;
        return frac_area;
    }
        // for 0 <= x < y:
    double r = sqrt(pow(x,2) + pow(y,2));
    double a = 0.5 + (x/2 + r*(r - R)/l)/y;
    if (a <= 0) // Case 1, pixel below circle
    {
        frac_area = 1;
    }
    else if (a >= 1 + x/y) // Case 5, pixel above circle
    {
        frac_area = 0;
    }
    else if ((x/y <= a) and (a <= 1)) // Case 3, circle intersects left 
                                      // and right sides
    {
        frac_area = 0.5 - (r - R)*r/(l*y);
    }
    else if ((0 < a) and (a < x/y)) // Case 2, circle intersects top 
                                    // and right
    {
        frac_area = 1.0 - 0.5*(y/x)*pow(a, 2);
    }
    else // Case 4: (1 < a) and (a < 1 + x/y), 
         // circle intersects left and bottom
    {
        frac_area = (x/(2*y))*pow((1 - (a - 1)*(y/x)), 2);
    }
    
    return frac_area;
}

double ill_frac(double x, double y, double Rin, double Rout, double l)
{
    // reflect x,y into pi/4 <= theta <= pi/2 region of pupil plane
    if(x<0)
        x *= -1;
    if(y<0)
        y *= -1;
    if(y < x)
    {
        const double tmp_holder = x;
        x = y;
        y = tmp_holder;
    }
    // initialize illumated fraction to case of pixel outside illuminated area
    double frac_area = 0.0;
    // compute fractional area if pixel inside or near border of illum area
    const double r = sqrt(pow(x,2) + pow(y,2));
        // pixel near Rin
    if((Rin - l/sqrt(2) < r) and (r < Rin + l/sqrt(2)))
        frac_area = 1 - frac_below_circle(x, y, Rin, l);
        // pixel near Rout
    if((Rout - l/sqrt(2) < r) and (r < Rout + l/sqrt(2)))
        frac_area = frac_below_circle(x, y, Rout, l);
        // pixel entirely in illuminated area
    if((Rin + l/sqrt(2) < r) and (r < Rout - l/sqrt(2)))
        frac_area = 1;
    return frac_area;
}