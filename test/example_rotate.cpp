#include <iostream>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

using namespace std;

int main( int argc, char* argv[] ){

    double alpha = 1E-3;
    gsl_vector_complex* evals = gsl_vector_complex_alloc(3);
    gsl_matrix_complex* evecs = gsl_matrix_complex_alloc(3,3);

    gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc(3);
    gsl_eigen_nonsymmv_params(1,w);
    
    double R_data[] = { cos(alpha), sin(alpha), 0.0,
                       -sin(alpha), cos(alpha), 0.0,
                        0.0       , 0.0       , 1.0 };

    gsl_matrix_view R = gsl_matrix_view_array(R_data, 3, 3);
    
    gsl_eigen_nonsymmv(&R.matrix, evals, evecs, w);

//  out put the results
    for( int i=0; i<3; ++i ){
        gsl_complex x = gsl_vector_complex_get( evals, i );
//        cout << "eval--> " << GSL_REAL(x) << "\t" << GSL_IMAG(x) << "\n";
//        cout << "norm--> " << sqrt(GSL_REAL(x)*GSL_REAL(x)*0 + GSL_IMAG(x)*GSL_IMAG(x)) << "\n";
        
        double x_real = GSL_REAL(x);
        double x_imag = GSL_IMAG(x);

//        if( (fabs(x_imag) < 1E-10) && ((fabs(x_real)-1.0) < 1E-10) ){

            cout << "eigen-values are: \n";
            cout << "(" << x_real << ", " << x_imag << ")" << endl;
            cout << "eigen-vector is:\n";
            for( int j=0; j<3; ++j ){
                gsl_complex y = gsl_matrix_complex_get( evecs, i, j );
                cout << "\tevec--> " << GSL_REAL(y) << "\t" << GSL_IMAG(y) << "\n";
            }

//        }
    }

    gsl_vector_complex_free(evals);
    gsl_matrix_complex_free(evecs);
    gsl_eigen_nonsymmv_free(w);
}
