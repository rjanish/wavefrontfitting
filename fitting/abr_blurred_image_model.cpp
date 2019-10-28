#include "abr_fit_image.h"
#include <cmath>
//DEBUG
// #include <iostream>
// #include <fstream>
// #include <cfloat>
//END DEBUG
using namespace std;

#define KERNEL_DIM 13
#define NUM_BLR_PARAMS 10
#define NUM_UNBLR_PARAMS 9

void blurred_image_model(double blr_x_im, double blr_y_im, 
                         double blr_params[],
                         int num_blr_params, 
                         double *blr_brightness, 
                         double blr_partials[], 
                         double pixel_width)
{
    // DEBUG
    // ofstream steps("param_sorting.txt");
    // steps.precision(DBL_DIG);
    // END DBUG


    // useful names
    const double seeing_rad = blr_params[5];
    const int num_unblr_params = num_blr_params - 1;
    double unblr_params[NUM_UNBLR_PARAMS + 1];
    for (int param = 1; param <= 4; param++)
    {
        unblr_params[param] = blr_params[param];
    }
    for (int param = 5; param <= num_unblr_params; param++)
    {
        unblr_params[param] = blr_params[param+1];
    }
    const int kcenter = KERNEL_DIM/2.0; // cast to int truncates 
                                        // (ie, floor function)

    //DEBUG
    // for (int param = 1; param <= num_blr_params; param++)
    // {
    //     steps << blr_params[param] << endl;
    // }
    // steps << endl;
    // for (int param = 1; param <= num_blr_params - 1; param++)
    // {
    //     steps << unblr_params[param] << endl;
    // }
    // steps << endl
    //       << seeing_rad << endl << endl;
    //END DEBUG
    

    // compute kernel and kernel coordinates
    double kernel[KERNEL_DIM][KERNEL_DIM];
    double kernel_sum = 0;
    int kernel_size = 0;
    for (int ky = 0; ky < KERNEL_DIM; ky++)
    {
        for (int kx = 0; kx < KERNEL_DIM; kx++)
        {
            double kr = sqrt(pow(kx - kcenter, 2) + \ 
                             pow(kcenter - ky, 2));
            double kernel_value = 1 - kr/seeing_rad;
            if (kernel_value < 0)
            {
                kernel_value = 0.0;
            }
            else
            {
                kernel_size += 1;
            }
            kernel[ky][kx] = kernel_value;
            kernel_sum += kernel_value;
        }    
    }
    // normalize kernel
    for (int ky = 0; ky < KERNEL_DIM; ky++)
    {
        for (int kx = 0; kx < KERNEL_DIM; kx++)
        {
            kernel[ky][kx] *= (1.0/kernel_sum);
        }    
    }

    // compute unblurred quantities brightness
    double unblr_brightness[KERNEL_DIM][KERNEL_DIM];
    double unblr_partials[KERNEL_DIM][KERNEL_DIM][NUM_UNBLR_PARAMS + 1];
    for (int ky = 0; ky < KERNEL_DIM; ky++)
    {
        for (int kx = 0; kx < KERNEL_DIM; kx++)
        {
            double x_unblr = blr_x_im + (kx - kcenter)*pixel_width;
            double y_unblr = blr_y_im - (ky - kcenter)*pixel_width;
            double unblr_brightness_tmp = 0;
            double unblr_partials_tmp[NUM_UNBLR_PARAMS + 1];
            unblurred_image_model(x_unblr, y_unblr, 
                                  unblr_params,
                                  num_unblr_params, 
                                  &unblr_brightness_tmp,
                                  unblr_partials_tmp, 
                                  pixel_width);
            unblr_brightness[ky][kx] = unblr_brightness_tmp;
            for (int param = 1; param <= num_unblr_params; param++)
            {
                unblr_partials[ky][kx][param] = \
                        unblr_partials_tmp[param];
            }
        }    
    }

    // compute blurred brightness and derivatives
    double blr_bightness_tmp = 0.0;
    double blr_partials_tmp[NUM_BLR_PARAMS + 1];
    for (int param = 0; param <= NUM_BLR_PARAMS; param++)
    {
        blr_partials_tmp[param] = 0.0;
    }
    for (int ky = 0; ky < KERNEL_DIM; ky++)
    {
        for (int kx = 0; kx < KERNEL_DIM; kx++)
        {
            blr_bightness_tmp += \
                kernel[ky][kx]*unblr_brightness[ky][kx];

            for (int param = 1; param <= num_unblr_params; param++)
            {
                blr_partials_tmp[param] += \
                    kernel[ky][kx]*unblr_partials[ky][kx][param];
            }

            double kernel_diff = 1 - kernel[ky][kx]*kernel_size;
            double deriv_kernel = \
                    (1.0/(kernel_sum*seeing_rad))*kernel_diff;
            blr_partials_tmp[num_unblr_params + 1] += \
                    deriv_kernel*unblr_brightness[ky][kx];
            

        }    
    }

    *blr_brightness = blr_bightness_tmp;
    for (int param = 1; param <= 4; param++)
    {
        blr_partials[param] = blr_partials_tmp[param];
    }
    blr_partials[5] = blr_partials_tmp[10];
    for (int param = 6; param <= 10; param++)
    {
        blr_partials[param] = blr_partials_tmp[param-1];
    }

    // DEBUG
    // steps.close();
    // END DBUG

    return;
}


#undef KERNEL_DIM
#undef NUM_BLR_PARAMS
#undef NUM_UNBLR_PARAMS
