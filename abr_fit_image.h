
#ifndef FIT_IMAGE_HEADER
#define FIT_IMAGE_HEADER

// functions   
void grad_wavefront_residual(int n, double x[], double f[], double params[], 
                             double x_im, double y_im);
    // Given coordinates on the defocus-only image x, a set of aberration parameters params, and a set of image plane coordinates x_im and y_im, this returns in f the residual vector between the given image point and the image point produced by the defocus-only point, calculated using the gradient of the wavefront. "n" is the dimension of x and f, which in this case is always 2. 
void unpack_fits_data(char* file_name, double flat_image[], double image_pos[], 
                      double image_neg[], int num_data_pts);
    // need to give filename the correct type
    // Read FITS file with name filename and stores both images into separate 1d image arrays, and stores a 1d 'flat_image' array that contains the data from both images concatenated. 
void compute_error(double flat_image[], double flat_errors[],
                   int num_data_pts);
    // Computes the typical statistical error in each pixel by computing the standard deviation in the non-illuminated pixels, assuming that their mean corresponds to zero intensity.
void lev_mar_fit(double pixel_number[], double flat_image[], double errors[],
                 int num_data_pts, double fit_params[], int num_fit_params,
                 double **covar_matrix, double *chisq, double chisq_tol, 
                 int num_successes, int verbose, int seeing);
    // This runs the Numerical Recipes Levenberg-Marquardt routine to fit the data in 'flat_image'. This array should contain the data from a pair of defocused images concatenated into a single 1d array. The data will be fit for the parameters in 'fit_params', which should contain the initial guesses when the function is called.  The first num_common_fit_params in the array fit_params will be fit simultaneously to both images, the next (num_common_fit_params - num_fit_params)/2 will be fit to only the positive defocus image and the final (num_common_fit_params - num_fit_params)/2 will be fit to only the negative defocus image. The fitting algorithm will iterate until the change in chi squared is less than chisq_tol for num_successes consecutive iterations. After this, the final values of the fit params, the covariance matrix and the final chi squared are returned. If verbose is nonzero, the result of each iteration will be printed to stdout. 
void flat_image_model(double pixel_number, double fit_params[], 
                      double *brightness, double derivatives[], 
                      int num_fit_params);
    // This is a wrapper for the function 'image_model' that allows the model to be accessed as though it were a one-variable function. Given a 'pixel_number', this function determines which of the two images this pixel is located in, and it determines it's coordinate in the image plane of that image. It then selects only those parameters from the params array that pertain to this image, and uses them to call 'image_model'. The derivatives with respect to the parameters that pertain to the other image are set to zero before the brightness and full set of derivatives are returned.  
void unblurred_image_model(double x_im, double y_im, double abr_params[],
                 int num_abr_params, double *brightness, double partials[],
                 double pixel_width);
    // Given an image coordinate and a set of aberration parameters, this computes the brightness and the partial derivatives with respect to the aberration parameters of the pixel centered on that coordinate. 
void blurred_image_model(double x_im, double y_im, double abr_params[],
                 int num_abr_params, double *brightness, double partials[],
                 double pixel_width);

int project_to_defocus_only(double x_im, double y_im, double params[],
                            double* x_df, double* y_df);
    // Using the given aberration parameters, this projects the image point x_im, y_im to the defocus-only image. This is done by finding the zero (i.e., the first zero located) of the 'grad_wavefont_residual' function using the Numerical Recipes 2d newton method routine. The return value will indicated whether or not the algorithm converged: 0 -> convergence, 1 -> no convergence. 
double inverse_hessian_det(double x_df, double y_df, double params[]);
    // returns the determinant of the inverse hessian of the wavefront, evaluated in the defocus-only image.
void inverse_hessian(double x_df, double y_df, double params[], 
                     double Hinv[][2], double Hinv_determinant);
    // returns the inverse of the Hessian matrix of the wavefront evaluated in the defocus plane.
void grad_df_position(double x_df, double y_df, double params[],
                      double Hinv[][2], double grad_x_df[], 
                      double grad_y_df[]); 
    // Assuming a fixed image plane point, this returns the gradient of the projected x and y defocus-only positions, with respect to the aberration parameters.    
void grad_bulk_term(double x_df, double y_df, double grad_x_df[],
                    double grad_y_df[], double Hinv[][2],
                    double params[], double bulk_term[]);
void grad_edge_term(double x_df, double y_df, double grad_x_df[],
                    double grad_y_df[], double params[], double l, 
                    double Rin, double Rout, double edge_term[]);
    // Compute the gradient of pixel brightness with respect to all aberration parameters. The computation is divided into two separate functions, the sum of which is proportional to the required gradient.  The edge term is only nonzero if the pixel lies on the edge of the defocus-only illuminated area. 
double ill_frac(double x_df, double y_df, double Rin,
                double Rout, double pixel_length);
    // Returns the fraction of a pixel with length pixel_length, and center x_df, y_df that is illuminated in the defocus-only image.
double frac_below_circle(double x_df, double y_df, 
                         double radius, double pixel_length);
    // Returns the fraction of the area of a pixel of length pixel_length and center x_df, y_df that is located inside the circle centered on the origin with the given radius. This uses the approximation that radius >> pixel_length, so that the circle can be replaced with its tangent line near the pixel.
double partial_frac_below_circle_wrt_defocus(double x, double y, double R,
                                             double l, double defocus);
    // Computes the partial derivative of the 'frac_below_circle' with respect to the defocus parameter. There is a direct dependence on the defocus because the radii of the illuminated area in the defocus-only plane depend on the defocus parameter. 
void grad_frac_below_circle_wrt_dfposition(double x, double y, double R,
                                           double frac_below_grad[2], 
                                           double l); 
    // Computes the gradient of the 'frac_below_circle' with respect to the defocus-only image coordinates.
void compute_initial_guesses(double flat_image[],
                             int num_fit_params, int num_common_fit_params,
                             int num_x_pixels, int num_y_pixels, 
                             double param_guesses[]);
    // Compute initial guesses for the parameters by analyzing each image for each aberration in turn, assuming that all other aberrations are zero. 

#endif 
