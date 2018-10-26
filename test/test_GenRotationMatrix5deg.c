#include "SurveySim.h"

int main( int argc, char* argv[] ){

    int p_rank;
    int p_size = 0;
    MPI_Init(NULL, NULL );
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p_size);

	if( argc < 5 ){
		printf("usage:\n   %s ra_old dec_old ra_new dec_new\n",argv[0]);
		exit(0);
	}

	double ra_old = atof(argv[1]);
	double dec_old= atof(argv[2]);
	double ra_new = atof(argv[3]);
	double dec_new= atof(argv[4]);

	double R5deg[9] = {0};
	double axis[3] = {0};
	int status = GenRotationMatrix5deg(ra_old,dec_old,ra_new,dec_new,R5deg,axis,p_rank);

	if( status == 0 ){
		printf("--- Successfuly generated the 5deg rotation matrix ---\n");
		print_mat("5Rdeg",R5deg);
	}
	else{
		print_Error_Msg("failed to generate the 5deg rotation matrix",p_rank);
	}

	MPI_Finalize();

	return 0;
}
