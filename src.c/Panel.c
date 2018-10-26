#include "SurveySim.h"

int TestPanelAngle(	double *sun, 
					double dist_sun,
					double x,
					double y,
					double z,
					double nx,
					double ny,
					double nz,
					double *cosval,
					double rotate_angle_max,
					int p_rank ){

	int status = 0;

	// return 1;

	//	望远镜指向
	double p[3] = {x,y,z};

	//	帆板法线初始指向
	double n_p[3] = {nx,ny,nz}; //已经归一化
	
	//	选择朝着太阳的那个法线
	//	(从与唐怀金的讨论中意识到有这种可能的情况;如果不进行法线的选择,那么计算得到的角度可能会大于90度)
	if( DotProduct(sun,n_p) < 0 ){
		n_p[0] = -nx;
		n_p[1] = -ny;
		n_p[2] = -nz;
	}

	//	帆板旋转轴:由望远镜指向与帆板法线初始指向叉乘得到
	double axis[3] = {0,0,0};
	CrossProduct(p,n_p,axis);	//无需再归一化


	//	开始对旋转角进行循环，搜索最小的角度，也就是最大的余弦值（其实是绝对值）
	double R[9];	// store elements of the rotation matrix
	double cosval_tmp = 0; // 保存余弦值，并逐步比较，找到最大值
	double n_p_tmp[3] = {0,0,0};	// 保存临时的旋转后的帆板法线

	double tmp = 0;
	double alpha_min = -1.0*fabs(rotate_angle_max);
	double alpha_max = fabs(rotate_angle_max);
	double alpha = alpha_min;
	
	while( alpha <= alpha_max ){
		
		GenRotationMatrix(axis,alpha,R,p_rank);
		MatrixVecProduct(R,n_p,n_p_tmp);
		tmp = DotProduct(sun,n_p_tmp)/dist_sun;

		if( fabs(tmp) > cosval_tmp ){
			cosval_tmp = fabs(tmp);
			cosval_tmp = (cosval_tmp > 1.0) ?  1.0 : cosval_tmp;
			*cosval = cosval_tmp;
		}
		
		if( fabs(cosval_tmp-1.0) < 1e-5 ){
			return 1;
		}

		alpha += 1;
	}

	if( cosval_tmp >= COS_SUN_PLANE_ANGLE ){
		status = 1;
	}

	return status;
}