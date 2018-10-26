//	Test transformations of coordinates from Ecliptic reference frame to Galactic reference frame.
//	All the transformation tools used in this main function are implemented by Zhang Xin

#include "SurveySim.h"

double gb_zhang( double ra, double dec){
	double bn = 57.2957795
				* asin(-0.8676660 * cos(ra * PI_180) * cos(dec * PI_180)
				- 0.1980764 * sin(ra * PI_180) * cos(dec * PI_180)
				+ 0.4559840 * sin(dec * PI_180) );
	return bn;
}

int main( int argc, char* argv[] ){

	if( argc < 3 ){
		printf("usage: %s L B\n",argv[0]);
		printf("note: (L,B) is coordinate in the Ecliptic reference frame, in unit of degrees\n");
		exit(0);
	}

	double L = atof(argv[1]);
	double B = atof(argv[2]);

	printf("Ecliptic (L,B) : %10.5f, %10.5f\n",L,B);

	L = L*M_PI/180;
	B = B*M_PI/180;

	double x_ec[3] = {cos(B)*cos(L), cos(B)*sin(L), sin(B)};	//	黄道坐标
	double x_eq[3];	//赤道坐标
	double x_gl[3];	//银道坐标

	//	从黄道坐标系变换到赤道坐标系
	// CoordinateSpinEquatorial2Ecliptic(x_eq,x_ec);
	CoordinateSpin_x(x_ec, x_eq, -23.4522);

//	for comparison to Zhang Xin's calculation
	double X_tmp[2];
	Cartesian2Equatorial(x_eq, X_tmp);
	double bn = gb_zhang(X_tmp[0],X_tmp[1]);

	//	从赤道坐标系到银道坐标系的旋转矩阵
	gsl_matrix *R_eq2gl = gsl_matrix_alloc(3,3);

	gsl_matrix_set(R_eq2gl,0,0,-0.05488);
	gsl_matrix_set(R_eq2gl,0,1,-0.87344);
	gsl_matrix_set(R_eq2gl,0,2,-0.48384);

	gsl_matrix_set(R_eq2gl,1,0,0.49411);
	gsl_matrix_set(R_eq2gl,1,1,-0.44483);
	gsl_matrix_set(R_eq2gl,1,2,0.74698);

	gsl_matrix_set(R_eq2gl,2,0,-0.86767);
	gsl_matrix_set(R_eq2gl,2,1,-0.19808);
	gsl_matrix_set(R_eq2gl,2,2,0.45598);


	//	从银道坐标系变换到赤道坐标系
	gsl_vector *x_old = gsl_vector_alloc(3);
	gsl_vector *x_new = gsl_vector_alloc(3);

	gsl_vector_set(x_old,0,x_eq[0]);
	gsl_vector_set(x_old,1,x_eq[1]);
	gsl_vector_set(x_old,2,x_eq[2]);

	gsl_blas_dgemv(CblasNoTrans,1,R_eq2gl,x_old,0,x_new);

	x_gl[0] = gsl_vector_get(x_new,0);
	x_gl[1] = gsl_vector_get(x_new,1);
	x_gl[2] = gsl_vector_get(x_new,2);


	//	从笛卡尔坐标变换为球坐标
	double X[2];
	Cartesian2Equatorial(x_gl, X);

	printf("Galactic (l,b) : %10.5f, %10.5f\n",X[0],X[1]);
	// printf("bn= %8.5f\n",bn);

	gsl_vector_free(x_old);
	gsl_vector_free(x_new);
	gsl_matrix_free(R_eq2gl);

}

