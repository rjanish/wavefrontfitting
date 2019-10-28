#include <cmath>
#include <cstring>
#include <iostream>
#include "nrutil.h"
#include "abr_fit_image.h"
using namespace std;

#define NUM_PIXELS_X 120
#define NUM_PIXELS_Y 120
#define TOTAL_PIXELS 2*120*120
#define NUM_FIT_PARAMS 15
#define NUM_COMMOM_FIT_PARAMS 5

int main(int argc, char* argv[])
{
    // Parameter Convention:
    // param array for single unblurred image = [empty, coma_x, coma_y,
    //                                           astig_x, astig_y, tilt_x,
    //                                           tilt_y, defocus,   
    //                                           pupil_flux, background]
    // total param array = [empty, coma_x, coma_y, astig_x, astig_y, 
    //                      tilt_x+, tilt_y+, defocus+, 
    //                      pupil_flux+, background+, seeing+,
    //                      tilt_x-, tilt_y-, defocus-, 
    //                      pupil_flux-, background-, seeing-] 
    // In both cases the paramters start at index 1 and index 0 holds 
    // garbage data. This is because the nr lev-mar routine wants to work 
    // with 1-indexed arrays. 

    // if -v argument is not used, the output is:
    //              0: FITS filename
    //        [1, 16]: initial guesses
    //       [17, 32]: best fit values
    //      [33, 288]: covariance matrix 
    //            289: number of iterations
    //            290: final chi squared value
    //            291: total degrees of freedom 
    //            292: final lamda value (lev-mar damping factor)

    int verbose = 0;
    int start_of_filenames = 1;
    if (strcmp(argv[1], "-v") == 0)
    {
        verbose = 1;
        start_of_filenames = 2;
    }

    for (int file_num = start_of_filenames; file_num < argc; file_num++)
    {
        if (verbose)
        {
            cout << "\nFitting image: " << argv[file_num] 
            << " ......." << endl;
        }
        else
        {
            cout << argv[file_num] << endl;
        }

        // read in fits data
        double flat_image_array[TOTAL_PIXELS];
        double image_pos[NUM_PIXELS_X*NUM_PIXELS_Y];
        double image_neg[NUM_PIXELS_X*NUM_PIXELS_Y];
        unpack_fits_data(argv[file_num], flat_image_array, image_pos, 
                         image_neg, TOTAL_PIXELS);

        // compute initial guesses
        // the nr levmar routine wants the inital guesses as an array
        // indexed from 1, so the 0 index of 'param_guesses' will hold 
        // garbage and the actual paramters are in index 1 to NUM_FIT_PARAMS.
        double param_guesses[NUM_FIT_PARAMS + 1]; 
        compute_initial_guesses(flat_image_array, NUM_FIT_PARAMS,
                                NUM_COMMOM_FIT_PARAMS, NUM_PIXELS_X, 
                                NUM_PIXELS_Y, param_guesses,
                                argv[file_num]);

        // compute errors
        double flat_errors_array[TOTAL_PIXELS];
        compute_error(flat_image_array, flat_errors_array, TOTAL_PIXELS);

        // make 'x values' vector and package data and error arrays into nr vectors
        double *pixel_number_vec, *flat_image_vec, *flat_errors_vec;
        pixel_number_vec = vector(1, TOTAL_PIXELS); 
        flat_image_vec = vector(1, TOTAL_PIXELS); 
        flat_errors_vec = vector(1, TOTAL_PIXELS); 
        for (int i = 1; i <= TOTAL_PIXELS; i++)
        {
            pixel_number_vec[i] = i-1;
            flat_image_vec[i] = flat_image_array[i-1];
            flat_errors_vec[i] = flat_errors_array[i-1];
        }

        // containers for fit results
        double **covar_matrix; 
        covar_matrix = matrix(1, NUM_FIT_PARAMS, 1, NUM_FIT_PARAMS);
        double chisq = 0;
        
        // fitting routine procedural parameters
        double chisq_tol = 0.01;
        int num_successes = 5;
        
        // set seeing guesses and re-fit
        int seeing = 0;
        lev_mar_fit(pixel_number_vec, flat_image_vec, flat_errors_vec,
                    TOTAL_PIXELS, param_guesses, NUM_FIT_PARAMS, 
                    covar_matrix, &chisq, chisq_tol, num_successes, 
                    verbose, seeing);
    
            
        
        // clean up
        free_vector(pixel_number_vec, 1, TOTAL_PIXELS); 
        free_vector(flat_image_vec, 1, TOTAL_PIXELS); 
        free_vector(flat_errors_vec, 1, TOTAL_PIXELS); 
        free_matrix(covar_matrix, 1, NUM_FIT_PARAMS, 1, NUM_FIT_PARAMS);
    

    }
}

#undef NUM_PIXELS_X
#undef NUM_PIXELS_Y
#undef TOTAL_PIXELS
#undef NUM_FIT_PARAMS
#undef NUM_COMMOM_FIT_PARAMS
