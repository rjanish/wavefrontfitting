#include "abr_fit_image.h"
#include "nr_routines.h"
#include "nrutil.h"
#include <cmath>
using namespace std;

#define NUM_PIXELS_X 120
#define NUM_PIXELS_Y 120
#define PIXEL_WIDTH 0.0000135
#define NUM_COMMON_FIT_PARAMS 4

#define DEFOCUS_SUM 8.525e-5

void flat_image_model(double pixel_number, double fit_params[], 
                      double *brightness, double gradient[], 
                      int num_fit_params)
{

    // Compute pixel coordinates 
        // The nr lev-mar routine wants the fitting function to have a 
        // double-type argument, but an integer type is needed to use the
        // mod function to compute the pixel coordinates
    int integer_pixel_number = pixel_number;
    int minus_defocus = integer_pixel_number/(NUM_PIXELS_X*NUM_PIXELS_Y);
    if (minus_defocus)
    {
        integer_pixel_number -= (NUM_PIXELS_X*NUM_PIXELS_Y);
    }
    int pixel_index_y = integer_pixel_number/NUM_PIXELS_X;
    int pixel_index_x = integer_pixel_number%NUM_PIXELS_X;
    double pixel_coordinate_x = \
                PIXEL_WIDTH*(-(NUM_PIXELS_X - 1.0)/2.0 + pixel_index_x);
    double pixel_coordinate_y = \
                PIXEL_WIDTH*(-(NUM_PIXELS_Y - 1.0)/2.0 + pixel_index_y);
    
    // Select applicable subset of parameters
    const int num_split_params = (num_fit_params - NUM_COMMON_FIT_PARAMS)/2;
    const int size_of_param_subset = NUM_COMMON_FIT_PARAMS + \
                                     num_split_params;
    double param_subset[size_of_param_subset + 1];
    for (int i = 1; i <= NUM_COMMON_FIT_PARAMS; i++)
    {
        param_subset[i] = fit_params[i];
    }
    const int start_of_unique_params = NUM_COMMON_FIT_PARAMS + \
                                       num_split_params*minus_defocus + 1;
    for (int j = start_of_unique_params; 
         j < (start_of_unique_params + num_split_params); 
         j++)
    {
        param_subset[j - num_split_params*minus_defocus] = fit_params[j];
    }
    

    // Eliminate independence of seeing and defocus terms, 
    // negate a_x, c_x, t_x for negative df image
    if (minus_defocus)
    {
        param_subset[7] = DEFOCUS_SUM - fit_params[7];
        param_subset[10] = fit_params[10];
        param_subset[1] = (-1)*param_subset[1]; // coma_x
        param_subset[3] = (-1)*param_subset[3]; // astig_x
        param_subset[5] = (-1)*param_subset[5]; // tilt_x
    }
    // call 2d image model 
    double gradient_subset[size_of_param_subset + 1];
    blurred_image_model(pixel_coordinate_x, pixel_coordinate_y, param_subset,
                       size_of_param_subset, brightness, gradient_subset,
                       PIXEL_WIDTH);
    
    // reconstruct full gradient
    for (int i = 1; i <= NUM_COMMON_FIT_PARAMS; i++)
    {
        gradient[i] = gradient_subset[i];
    }
    const int distance_to_dual_param = num_split_params*(1-2*minus_defocus);
    for (int j = start_of_unique_params; 
         j < (start_of_unique_params + num_split_params); 
         j++)
    {
        gradient[j] = gradient_subset[j - num_split_params*minus_defocus];
        gradient[j + distance_to_dual_param] = 0.0;
    }


    return;
}

#undef NUM_PIXELS_X
#undef NUM_PIXELS_Y
#undef PIXEL_WIDTH 
#undef NUM_COMMOM_FIT_PARAMS
#undef DEFOCUS_SUM