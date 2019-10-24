#include <cmath>
#include <iostream>
#include "abr_fit_image.h"
#include "nr_routines.h"
#include "nrutil.h"
using namespace std;

void lev_mar_fit(double pixel_number[], double flat_image[], double errors[],
                 int num_data_pts, double fit_params[], int num_fit_params,
                 double **covar_matrix, double *chisq, 
                 double chisq_tol, int num_successes, 
                 int verbose, int seeing)
{
    cout.precision(15);
    // print guesses
    if (verbose and seeing)
    {
        cout << "Inital guesses:" << endl;
        for (int i = 1; i <= num_fit_params; i++) 
            {
                cout << "\tfit_params[" << i << "] = " 
                     << fit_params[i] << endl;
            }
    }
    else
    {
        for (int i = 1; i <= num_fit_params; i++) 
            {
                cout << fit_params[i] << endl;
            }
    }

    // initalize
    double lamda = -1; // -1 signals nr routine to initialize
    double current_chisq = 0; 

    int *param_selector; // indicates which parameters are to be fitted for 
    param_selector = ivector(1, num_fit_params);
    
    for (int i = 1; i <= num_fit_params; i++)
    {
        param_selector[i] = 1;
    }
    param_selector[10] = seeing;
    param_selector[16] = seeing;

    double **truncated_hessian_matrix; // working space for nr mrqmin
    truncated_hessian_matrix = matrix(1, num_fit_params, 1, num_fit_params);

    mrqmin(pixel_number, flat_image, errors, num_data_pts, fit_params,
           param_selector, num_fit_params, covar_matrix, 
           truncated_hessian_matrix, &current_chisq, 
           flat_image_model, &lamda, verbose); // initializes and does 1 iteration
    
    // iterate fitting procedure
    double old_chisq = 0;
    double dof = num_data_pts - num_fit_params;
    int consecutive_successes = 0; 
    int iteration_number = 1;
    while (consecutive_successes < num_successes)
    {   

        // print results of previous iteration
        if (verbose and seeing)
        {
            cout << "Iteration " << iteration_number << ":" << endl;
            cout << "\tnew chi squared = " << current_chisq << endl;
            cout << "\tnew chisq/dof = " << current_chisq/dof << endl;
            if (iteration_number > 1)
            {
                cout << "\tdelta(chi squared) = " 
                        << current_chisq - old_chisq << endl;
                cout << "\tdelta(chisq/dof) = " 
                        << (current_chisq - old_chisq)/dof << endl;
            }
            cout << "\tnew lamda = " << lamda << endl;
            for (int i = 1; i <= num_fit_params; i++) 
            {
                cout << "\tfit_params[" << i << "] = " 
                     << fit_params[i] << endl;
            }
        }

        // iterate once 
        old_chisq = current_chisq;
        mrqmin(pixel_number, flat_image, errors, num_data_pts, fit_params,
               param_selector, num_fit_params, covar_matrix, 
               truncated_hessian_matrix, &current_chisq, 
               flat_image_model, &lamda, verbose);        
        
        // check for convergence and adjust iterration count
        if (current_chisq > old_chisq)
        {
            consecutive_successes = 0;
        }
        else if (abs(old_chisq - current_chisq) < 0.01)
        {
            consecutive_successes++;
        }
        iteration_number++;
    }

    // get final values
    double clean_up_lamda = 0; // lambda = 0 initiates clean-up
    mrqmin(pixel_number, flat_image, errors, num_data_pts, fit_params,
           param_selector, num_fit_params, covar_matrix, 
           truncated_hessian_matrix, &current_chisq, 
           flat_image_model, &clean_up_lamda, verbose);     
    *chisq = current_chisq;
    // print final resutls
    if (verbose and seeing)
    {
        cout << " ********** final values ********** " << endl;
        cout << "Iterations: " << iteration_number << endl;
        cout << "\tfinal chi squared = " << current_chisq << endl;
        cout << "\tfinal chisq/dof = " << current_chisq/dof << endl;
    
        cout << "\tdelta(chi squared) = " 
                        << current_chisq - old_chisq << endl;
         cout << "\tdelta(chisq/dof) = " 
                        << (current_chisq - old_chisq)/dof << endl;
        cout << "\t final lamda = " << lamda << endl;
        cout << "\t covariance matrix = " << endl;
        for (int i = 1; i <= num_fit_params; i++)
        {
            cout << "\t\t";
            for (int j = 1; j <= num_fit_params; j++)
            {
                cout << covar_matrix[i][j] << "  ";
            }
            cout << endl;
        }
        for (int i = 1; i <= num_fit_params; i++) 
        {
            cout << "fit_params[" << i << "] = " 
                 << fit_params[i] 
                 << " +- " << sqrt(covar_matrix[i][i]) << endl;
        }
    }
    else if (seeing)
    {
        // best fit results
        for (int i = 1; i <= num_fit_params; i++) 
        {
            cout << fit_params[i] << endl;
        }
        // lev-mar covariance matrix
        for (int i = 1; i <= num_fit_params; i++)
        {
            for (int j = 1; j <= num_fit_params; j++)
            {
                cout << covar_matrix[i][j] << endl;
            }
        }
    }
    // fit convergence stats
    cout << iteration_number << endl
         << current_chisq << endl
         << dof << endl
         << lamda << endl;
    

    // clean up
    free_ivector(param_selector, 1, num_fit_params);
    free_matrix(truncated_hessian_matrix, 1, num_fit_params, 
                                          1, num_fit_params);
    return;
}
