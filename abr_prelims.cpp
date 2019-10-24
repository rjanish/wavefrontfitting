#include "abr_fit_image.h"
#include "fitsio.h"
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

#define NUM_PIXELS_Y 120
#define NUM_PIXELS_X 120
#define BITS_PER_PIXEL 16
#define PIXEL_WIDTH 0.0000135
#define R_IN 0.8251
#define R_OUT 1.85

void unpack_fits_data(char* file_name, double flat_image[], 
                      double image_pos[], double image_neg[], int num_data_pts)
{
    // Open FITS file and verify correct image format
    fitsfile *fptr;   // FITS file pointer, defined in fitsio.h from CFITSIO
    int status = 0;   // CFITSIO status value must be initialized to zero
    fits_open_file(&fptr, file_name, READONLY, &status);
    int num_hdus = 0;
    fits_get_num_hdus(fptr, &num_hdus, &status);
    if (num_hdus != 2)
    {
        cout << "Error: FITS must contain exactly two images" << endl;
        return;
    }
    for (int ext = 1; ext < 3; ext++)
    {
        int max_dims = 2;
        int type = 3; // 0 -> image, 1 -> ascii table, 2 -> binary table
        int bitpix = 0, naxis = 0;
        long naxes[2] = {0, 0};
        fits_movabs_hdu(fptr, ext, &type, &status);
        fits_get_img_param(fptr, max_dims, &bitpix, &naxis, naxes, &status);
        int correct_image = (type == 0)*(bitpix == BITS_PER_PIXEL)*(naxis == 2)
                           *(naxes[0] == NUM_PIXELS_Y)*(naxes[1] == NUM_PIXELS_X);
        if (!correct_image or status)
        {
            cout << "Error: Invalid image format" << endl
                 << "\ttype: " << type << endl
                 << "\tbitpix: " << bitpix << endl
                 << "\tnaxis: " << naxis << endl
                 << "\tnaxes[0]: " << naxes[0] << endl
                 << "\tnaxes[1]: " << naxes[1] << endl
                 << "\tstatus: " << status << endl;
                 
            return;
        }
    }
    // get data
    long starting_pixel[2] = {1, 1};
    int any_nul_values = 0;
    // save positive defocus single array
    fits_movabs_hdu(fptr, 1, NULL, &status);
    fits_read_pix(fptr, TDOUBLE, starting_pixel, NUM_PIXELS_X*NUM_PIXELS_Y, 
                  NULL, image_pos, &any_nul_values, &status);
    // save negative defocus single array
    fits_movabs_hdu(fptr, 2, NULL, &status);
    fits_read_pix(fptr, TDOUBLE, starting_pixel, NUM_PIXELS_X*NUM_PIXELS_Y, 
                  NULL, image_neg, &any_nul_values, &status);
    // close FITS
    fits_close_file(fptr, &status);
    // report any CFITSIO errors
    if (status)
    {
        char cfitsio_error_txt[30];
        fits_get_errstatus(status, cfitsio_error_txt);
        cout << cfitsio_error_txt << endl;
        return;
    }
    // make combined flat_image array
    for (int pixel = 0; pixel < num_data_pts; pixel++)
    {
        if (pixel < num_data_pts/2)
        {
            flat_image[pixel] = image_pos[pixel];
        }
        else 
        {
            flat_image[pixel] = image_neg[pixel - num_data_pts/2];
        }
    }
    return;
}
    
void compute_error(double flat_image[], double flat_errors[], 
                   int num_data_pts)
{
    for (int i = 0; i < num_data_pts; i++)
    {
        flat_errors[i] = 2*sqrt(flat_image[i]);
    }
    return;
}

void compute_initial_guesses(double flat_image[],
                             int num_fit_params, int num_common_fit_params,
                             int num_x_pixels, int num_y_pixels, 
                             double param_guesses[])
{
    const double PI = atan(1.0)*4;
    // gaussian mask width reaches 3/4 distance to image edge
    double msk_width[2] = \
        {5.0*PIXEL_WIDTH*(NUM_PIXELS_X+NUM_PIXELS_Y)/(32.0),
         5.0*PIXEL_WIDTH*(NUM_PIXELS_X+NUM_PIXELS_Y)/(32.0)}; 
    double bkgnd[2] = {0.0, 0.0};
    double tilt_x[2] = {0.0, 0.0};
    double tilt_y[2] = {0.0, 0.0};
    double defocus[2] = {0.0, 0.0};
    double flux[2] = {0.0, 0.0};
    for (int df_label = 0; df_label < 2; df_label++)
    {
        int iteration = 0;
        double delta_tilt = 0.0;
        double delta_msk_width = 0.0;
        const int max_iters = 20;
        do 
        {   double total_msk = 0.0;
            double total_i_msk = 0.0;
            double total_i = 0.0;
            double total_x_i_msk = 0.0;
            double total_y_i_msk = 0.0;
            double total_rsq_msk = 0.0;
            for (int pixel_num = 0; 
                 pixel_num < NUM_PIXELS_Y*NUM_PIXELS_X; 
                 pixel_num++)
            {
                // get cartesian coordinates
                const int y_pix = pixel_num/NUM_PIXELS_X;
                const int x_pix = pixel_num%NUM_PIXELS_X;
                const double x = \
                    PIXEL_WIDTH*(-(NUM_PIXELS_X - 1.0)/2.0 + x_pix);
                const double y = \
                    PIXEL_WIDTH*(-(NUM_PIXELS_Y - 1.0)/2.0 + y_pix);
                const double rsq = pow(x-tilt_x[df_label], 2) + \
                                   pow(y-tilt_y[df_label], 2);
                // compute mask and star intensity
                const double mask = exp(-rsq/pow(msk_width[df_label], 2));
                const int flat_index = pixel_num + \
                                 NUM_PIXELS_X*NUM_PIXELS_Y*df_label;
                const double i_star = \
                    flat_image[flat_index] - bkgnd[df_label];
                // compute moments
                total_msk = total_msk + mask;
                total_i_msk = total_i_msk + i_star*mask;
                total_i = total_i + i_star;
                total_x_i_msk = total_x_i_msk + i_star*x*mask;
                total_y_i_msk = total_y_i_msk + i_star*y*mask;
                total_rsq_msk = total_rsq_msk + i_star*rsq*mask;
            }
            // new tilt
            const double new_tilt_x = total_x_i_msk/total_i_msk;
            const double new_tilt_y = total_y_i_msk/total_i_msk;
            const double delta_tilt_x = new_tilt_x - tilt_x[df_label];
            const double delta_tilt_y = new_tilt_y - tilt_y[df_label];
            delta_tilt = sqrt(pow(delta_tilt_x, 2) + \
                                           pow(delta_tilt_y, 2));
            tilt_x[df_label] = new_tilt_x;
            tilt_y[df_label] = new_tilt_y;
            // new mask width
            const double avg_rad_sq = total_rsq_msk/total_i_msk;
            const double new_msk_width = sqrt(avg_rad_sq);
            delta_msk_width = abs(new_msk_width - msk_width[df_label]);
            msk_width[df_label] = new_msk_width;
            // new defocus, background and star flux
            defocus[df_label] = \
                sqrt(avg_rad_sq/(2*(pow(R_OUT, 2) + pow(R_IN, 2))));
            flux[df_label] = \
                total_i*pow(PIXEL_WIDTH,2)/(PI*(pow(R_OUT,2) - pow(R_IN, 2)));
            bkgnd[df_label] += \
                (total_i - total_i_msk)/(NUM_PIXELS_X*NUM_PIXELS_Y - total_msk);

            // check for infinite looping
            iteration = iteration + 1;
            if (iteration > max_iters)
            {
                break;
            }
        } while ((delta_tilt > PIXEL_WIDTH) and \
                 (delta_msk_width > PIXEL_WIDTH));
    }

    param_guesses[0] = 0.0;                 // garbage
    param_guesses[1] = 0.0;                 // coma x
    param_guesses[2] = 0.0;                 // coma y
    param_guesses[3] = 0.0;                 // asitg x
    param_guesses[4] = 0.0;                 // astig y
    param_guesses[5] = tilt_x[0];           // +
    param_guesses[6] = tilt_y[0];           // +  
    param_guesses[7] = defocus[0];          // +
    param_guesses[8] = flux[0];             // +
    param_guesses[9] = bkgnd[0];            // +
    param_guesses[10] = 0.5;                // + seeing
    param_guesses[11] = tilt_x[1];          // -
    param_guesses[12] = tilt_y[1];          // -
    param_guesses[13] = defocus[0];         // - this will become (d+ - C)
    param_guesses[14] = flux[1];            // -
    param_guesses[15] = bkgnd[1];           // -
    param_guesses[16] = 0.5;                // - this will become seeing+

    return;
}

#undef NUM_PIXELS_Y
#undef NUM_PIXELS_X 
#undef BITS_PER_PIXEL
#undef PIXEL_WIDTH
#undef R_OUT
#undef R_IN
