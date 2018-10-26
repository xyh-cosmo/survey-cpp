//	Test transformations of coordinates from Ecliptic reference frame to Galactic reference frame.
//	All the transformation tools used in this main function are implemented by Zhang Xin

#include "SurveySim.h"


int GetMatrixInverse( gsl_matrix* m, gsl_matrix* inv_m ){
	int i,j,signum;
	gsl_matrix *LU = gsl_matrix_alloc(3,3);
	gsl_permutation *p = gsl_permutation_alloc(3);

	for( i=0; i<3; i++ ){
		for( j=0; j<3; j++){
			gsl_matrix_set(LU,i,j,gsl_matrix_get(m,i,j));
		}
	}

	gsl_linalg_LU_decomp(LU,p,&signum);
	gsl_linalg_LU_invert(LU,p,inv_m);

	gsl_matrix_free(LU);
	gsl_permutation_free(p);

	return 0;
}

int main( int argc, char* argv[] ){

	if( argc < 3 ){
		printf("usage: %s l b\n",argv[0]);
		printf("note: (l,b) is coordinate in the Galactic reference frame, in unit of degrees\n");
		exit(0);
	}

	double l = atof(argv[1]);
	double b = atof(argv[2]);

	printf("Galactic (l,b) : %10.5f, %10.5f\n",l,b);

	l = l*M_PI/180;
	b = b*M_PI/180;

	double x_gl[3] = {cos(b)*cos(l), cos(b)*sin(l), sin(b)};	//	银道坐标
	// printf("x_gl = %8.5f, y_gl = %8.5f, z_gl = %8.5f\n",x_gl[0],x_gl[1],x_gl[2]);
	
	double x_eq[3];	//赤道坐标
	double x_ec[3];	//黄道坐标

	//	从赤道坐标系到银道坐标系的旋转矩阵
	gsl_matrix *R_eq2gl = gsl_matrix_alloc(3,3);
	//	从银道坐标系到赤道坐标系的旋转矩阵
	gsl_matrix *R_gl2eq = gsl_matrix_alloc(3,3);

	gsl_matrix_set(R_eq2gl,0,0,-0.05488);
	gsl_matrix_set(R_eq2gl,0,1,-0.87344);
	gsl_matrix_set(R_eq2gl,0,2,-0.48384);

	gsl_matrix_set(R_eq2gl,1,0,0.49411);
	gsl_matrix_set(R_eq2gl,1,1,-0.44483);
	gsl_matrix_set(R_eq2gl,1,2,0.74698);

	gsl_matrix_set(R_eq2gl,2,0,-0.86767);
	gsl_matrix_set(R_eq2gl,2,1,-0.19808);
	gsl_matrix_set(R_eq2gl,2,2,0.45598);

	GetMatrixInverse(R_eq2gl,R_gl2eq);

	//	从银道坐标系变换到赤道坐标系
	gsl_vector *x_old = gsl_vector_alloc(3);
	gsl_vector *x_new = gsl_vector_alloc(3);

	gsl_vector_set(x_old,0,x_gl[0]);
	gsl_vector_set(x_old,1,x_gl[1]);
	gsl_vector_set(x_old,2,x_gl[2]);

	gsl_blas_dgemv(CblasNoTrans,1,R_gl2eq,x_old,0,x_new);

	x_eq[0] = gsl_vector_get(x_new,0);
	x_eq[1] = gsl_vector_get(x_new,1);
	x_eq[2] = gsl_vector_get(x_new,2);

	//	从赤道坐标系变换到黄道坐标系
	CoordinateSpinEquatorial2Ecliptic(x_eq,x_ec);

	// printf("x_eq[0] = %8.5f, x_eq[1] = %8.5f, x_eq[2] = %8.5f\n",x_eq[0],x_eq[1],x_eq[2]);

	//	从笛卡尔坐标变换为球坐标
	double X[2];
	Cartesian2Equatorial(x_ec, X);

	printf("Ecliptic (L,B) : %10.5f, %10.5f\n",X[0],X[1]);

	gsl_vector_free(x_old);
	gsl_vector_free(x_new);
	gsl_matrix_free(R_eq2gl);
	gsl_matrix_free(R_gl2eq);

}

