#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include "nr_newt.h"
#include "nrutil.h"
using namespace std;

#define KERNEL_DIM 13

// telescope parameters
double rIn, rOut, pixel_width, flux;
int num_x, num_y;
// gradient parameters
double gradx_constant, gradx_x_lin, gradx_y_lin;
double gradx_x_sq, gradx_y_sq, gradx_cross;
double grady_constant, grady_x_lin, grady_y_lin;
double grady_x_sq, grady_y_sq, grady_cross;
// image pixel coordinates
double px, py;
// overall tilt + field coords displacment
double xShift, yShift;
// hessian determinatant parameters
double hd_constant, hd_x_lin, hd_y_lin, hd_x_sq, hd_y_sq, hd_cross;
// defoucs
double d;

// illuminated fraction functions
double frac_below_circle(double x, double y, double R, double l)
{
    /* For a circle of radius R centered on the orgin and a square of side
    length l and center at cartesian coords x and y, this function returns
    the fraction of the area of the square that is interior to the circle.

    It is assumed that R >> l, in order that the interior area may be computed
    by approximating the circle with its tangent line. The tangent line used is
    that wich intersects the circle closest to the point x, y.

    The coords x and y must satisify 0 <= x <= y (i.e., the center of
    the square is located at a polar angle of pi/4 <= theta <= pi/2) for
    this function to return the correct value. If this condition is not met,
    the square should be reflected about the x-axis and/or y-axis and/or line
    y=x as needed before the function is called. Reflecting the square in this
    way does not change the magnitude of the area interior to the circle. */

    double frac_area = 0;
        // for y = x = 0:
    if(y == 0)
    {
        frac_area = 1.0;
        return frac_area;
    }
        // for 0 <= x < y:
    double r = sqrt(pow(x,2) + pow(y,2));
    double a = 0.5 + x/(2*y) + ((r - R)/l)*(r/y);
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
        frac_area = 1.0 + x/(2*y) - a;
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
        double tmp_holder = x;
        x = y;
        y = tmp_holder;
    }
    // initialize illumated fraction to case of pixel outside illuminated area
    double frac_area=0;
    // compute fractional area if pixel inside or near border of illum area
    double r = sqrt(pow(x,2) + pow(y,2));
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

// grad(W) transformation system of equations
void funcv(int n, double x[], double f[])
{
    f[1] = (gradx_x_sq*x[1] + gradx_x_lin)*x[1] + \
           (gradx_y_sq*x[2] + gradx_y_lin + gradx_cross*x[1])*x[2] - px;

    f[2] = (grady_x_sq*x[1] + grady_x_lin)*x[1] + \
           (grady_y_sq*x[2] + grady_y_lin + grady_cross*x[1])*x[2] - py;
}

double inverse_hessian_determinant(double x, double y)
{
    const double hd = (hd_x_sq*x + hd_x_lin)*x + \
                      (hd_y_sq*y + hd_y_lin + hd_cross*x)*y + \
                       hd_constant;
    return (1.0/abs(hd));
}

void seeing_kernel(double r_0, int d, double k[][KERNEL_DIM])
{
    // compute convolution kernel for seeing
    double distance_sum = 0;
    double k_sum = 0;
    const int center = d/2.0; // cast to int -> floor
    for (int y = 0; y < d; y++)
    {
        for (int x = 0; x < d; x++)
        {
           double r = sqrt(pow((y - center), 2) + pow((x - center), 2));
           distance_sum += r;
           double kernel_value = 1 - r/r_0;
           if (kernel_value < 0)
           {
                kernel_value = 0.0;
           }
           k[y][x] = kernel_value;
           k_sum += kernel_value;
        }    
    }
    // normalize
    for (int y = 0; y < d; y++)
    {
        for (int x = 0; x < d; x++)
        {
            k[y][x] *= (1.0/k_sum);
        }    
    }
    return;
}

// N = number of equations & number of unknowns in transformation
#define N 2

// driver
int main(int argc, char* argv[])
{
    // read in telescope parameters
    cin >> num_x;
    cin >> num_y;
    double telescope_rad_in, telescope_rad_out;
    cin >> telescope_rad_in;
    cin >> telescope_rad_out;
    cin >> pixel_width;
    cin >> flux;
    // read in aberration parameters
    double coma, phi_coma, astig, phi_astig, defocus;
    double tilt, phi_tilt, sigma, theta, background, seeing;
    cin >> coma;
    cin >> phi_coma;
    cin >> astig;
    cin >> phi_astig;
    cin >> defocus;
    cin >> tilt;
    cin >> phi_tilt;
    cin >> sigma;
    cin >> theta;
    cin >> background;
    cin >> seeing;
    // set defocus plane radii
    rIn = abs(telescope_rad_in*2.0*defocus);
    rOut = abs(telescope_rad_out*2.0*defocus);
    // store defoucs in global
    d = defocus;
    // tilt and field coordinate displacment
    xShift = tilt*cos(phi_tilt) + sigma*cos(theta);
    yShift = tilt*sin(phi_tilt) + sigma*sin(theta);
    // compute global paramters for gradient transform
        // x-component
    gradx_x_lin = 1.0 + astig*cos(phi_astig)/defocus;
    gradx_y_lin = astig*sin(phi_astig)/defocus;
    gradx_x_sq = 3*coma*cos(phi_coma)/(4*pow(defocus, 2));
    gradx_y_sq = coma*cos(phi_coma)/(4*pow(defocus, 2));
    gradx_cross = coma*sin(phi_coma)/(2*pow(defocus, 2));
        // y-component
    grady_x_lin = astig*sin(phi_astig)/defocus;
    grady_y_lin = 1.0 - astig*cos(phi_astig)/defocus;
    grady_x_sq = coma*sin(phi_coma)/(4*pow(defocus, 2));
    grady_y_sq = 3*coma*sin(phi_coma)/(4*pow(defocus, 2));
    grady_cross = coma*cos(phi_coma)/(2*pow(defocus, 2));
    // compute global paramters for hessian determinant
    hd_constant = 1 - pow(astig, 2)/pow(defocus, 2);
    hd_x_lin = (coma/pow(defocus, 2))*(2*cos(phi_coma) - \
               (astig/defocus)*cos(phi_coma - phi_astig));
    hd_y_lin = (coma/pow(defocus, 2))*(2*sin(phi_coma) + \
               (astig/defocus)*sin(phi_coma - phi_astig));
    hd_x_sq = (pow(coma, 2)/(4*pow(defocus, 4)))*(1 + 2*cos(2*phi_coma));
    hd_y_sq = (pow(coma, 2)/(4*pow(defocus, 4)))*(1 - 2*cos(2*phi_coma));
    hd_cross = (pow(coma, 2)/pow(defocus, 4))*sin(2*phi_coma);

    // set output format to max precision
    cout.precision(DBL_DIG);
    // guess point for projection to pupil plane
    double *x = vector(1,N);
    // vector to hold im-to-pupil transfomation functions
    double *f = vector(1,N);
    // flags to indicate nonconvergence (check) and
    // convergence location (where)
    int check, where=0;


    // seeing kernel
    double kernel[KERNEL_DIM][KERNEL_DIM];
    seeing_kernel(seeing, KERNEL_DIM, kernel);

    // main loop
    for(double yIm = pixel_width*(num_y-1.0)/2.0 - yShift; yIm > -pixel_width*num_y/2.0 - yShift; yIm -= pixel_width)
    {
        for(double xIm = -pixel_width*(num_x-1.0)/2.0 - xShift; xIm < pixel_width*num_x/2.0 - xShift; xIm += pixel_width)
        {
            double brightness = 0;
            double old_brightness = 0;
            double kernel_sum = 0;

            for (int ky = 0; ky < KERNEL_DIM; ky++)
            {
                for (int kx = 0; kx < KERNEL_DIM; kx++)
                {
                    int base = KERNEL_DIM/2.0;
                    int shiftx = (-1)*base + kx;
                    int shifty = base - ky;
                    px = xIm + shiftx*pixel_width;
                    py = yIm + shifty*pixel_width;
                    double unblured = 0;

                    // set guess to image coords
                    x[1] = px;
                    x[2] = py;
                    // project back to pupil plane
                    newt(x, N, &check, &where, funcv);
                    // store pupil coords in x
                    funcv(N, x, f);

                    if(check) // no pupil point
                    {
                        // set brightness to background
                        unblured = background;
                    }
                    else // found a pupil point
                    {
                        // compute brightness
                        double defocus_plane_area = inverse_hessian_determinant(x[1], x[2]);
                        double illuminated_fraction = ill_frac(x[1], x[2], rIn, rOut, pixel_width);
                        double defocus_plane_flux = flux/(4.0*pow(d, 2));
                        
                        unblured = defocus_plane_flux*defocus_plane_area*illuminated_fraction + background;
                    }

                    brightness += unblured*kernel[ky][kx];
                    
                }               
            }

            // write brightness to stdout
            cout << brightness << " ";
        }
        // end current row
        cout << endl;

    }

    // release guess and transformation vectors from memory
    free_vector(f,1,N);
    free_vector(x,1,N);
    return 0;
}

#undef KERNEL_DIM
#undef N
