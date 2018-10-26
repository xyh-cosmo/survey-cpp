 // TODO: 修改计算转动角的那些函数的接口

 #include "SurveySim.h"

void ErrorInfo( char* info ){
    // printf("FILE: %s\n",__FILE__);
    // printf("LINE: %d\n",__LINE__);
    // printf("FunName: %s\n",__FUNCTION__);
    printf("\n**** ErrInfo: %s ****\n\n",info);
}

void print_vec(char* vecname, double* vec){
    printf("Vector Name: %s \n\t[",vecname);
    int i;
    for(i=0;i<3;i++)
        printf(" %8.5f ",vec[i]);
    printf("]\n\n");
}

void print_vec_diff(char* vecname1, double* vec1, char* vecname2, double* vec2){
    printf("difference between Vector '%s' and '%s' \n\t[",vecname1,vecname2);
    int i;
    for(i=0;i<3;i++)
        printf(" %8.5f ",vec2[i]-vec1[i]);
    printf("]\n\n");
}

void print_mat(char* matname, double* mat){
    printf("Matrix Name: %s\n",matname);
    int i,j;
    for(i=0;i<3;i++){
        printf("\t|");
        for(j=0;j<3;j++)
            printf(" %8.5f ",mat[3*i+j]);
        printf("|\n");
    }
    printf("\n");
}

void print_rotation_axis(char* matname, double* mat){
    double axis[3] = {0,0,0};
    axis[0] = mat[7]-mat[5];
    axis[1] = mat[2]-mat[6];
    axis[2] = mat[3]-mat[1];
    printf("rotation axis of : %s\n\t[",matname);
    int i;
    for(i=0;i<3;i++)
        printf(" %8.5f ",axis[i]);
    printf("]\n\n");
}

void get_rotation_axis(double* mat, double *axis){
    axis[0] = mat[7]-mat[5];
    axis[1] = mat[2]-mat[6];
    axis[2] = mat[3]-mat[1];
}

void print_rotation_angle(char* matname, double* mat){
    double angle = GetRotationAngleFromMatrix(mat);
    printf("rotation angle of %s is %g (deg)\n\n",matname,angle);
}

//  计算两个矢量的点乘（内积）
double DotProduct( double *u, double *v ){
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}

//  计算两个矢量的叉乘
int CrossProduct( double *u, double *v, double *w ){

    w[0] = u[1]*v[2]-u[2]*v[1];
    w[1] =-u[0]*v[2]+u[2]*v[0];
    w[2] = u[0]*v[1]-u[1]*v[0];

    return 0;
}

double VecNorm( double *v ){
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double VecDiffNorm( double *u, double *v ){
	return sqrt( (u[0]-v[0])*(u[0]-v[0]) + (u[1]-v[1])*(u[1]-v[1]) + (u[2]-v[2])*(u[2]-v[2]) );
}

double GetMatrixTrace( double *M ){
    double trace = 0;
    int i=0;
    for( i=0; i<3; i++ ){
        trace += M[3*i+i];
    }

    return trace;
}

double GetRotationAngleFromMatrix(double *M){
    double trace = GetMatrixTrace(M);
    double tmp = 0.5*(trace-1);
    tmp = (tmp > 1.0) ? 1.0 : tmp;
    tmp = (tmp <-1.0) ?-1.0 : tmp;
    double angle = 180.0/PI*acos(tmp);

    if( gsl_isnan(angle) )
        angle = 180;
    return angle;
}

//  计算矩阵与向量的乘积: W=M*V
int MatrixVecProduct( double *M, double *V, double *W ){
    int i,j;
    for( i=0; i<3; i++ ){
        W[i] = 0;
        for( j=0; j<3; j++ ){
            W[i] += M[3*i+j]*V[j];
        }
    }

    return 0;
}

//  计算矩阵的乘积：T = R*S
int MatrixMultiplication( double* R, double* S, double* T ){

    // printf("==== debug inside MatrixMultiplication ====\n");
    // print_mat("R",R);
    // print_mat("S",S);

    int i,j,k;
    for( i=0; i<3; i++ ){
        for( j=0; j<3; j++ ){
            T[3*i+j] = 0;
            for( k=0; k<3; k++ ){
                T[3*i+j] += R[3*i+k]*S[3*k+j];
            }

            // printf("T[%d,%d] = %g\n",i,j,T[3*i+j]);
        }
    }

    return 0;
}

//  根据指定的转动轴和转动角度生成旋转矩阵（这里的旋转矩阵保存在1维的指针中）
int GenRotationMatrix( double* u, double angle_deg, double* R, int p_rank ){
    double theta = angle_deg*PI_180;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double One_cos_theta = 1-cos_theta;
    double ux = u[0];
    double uy = u[1];
    double uz = u[2];

    double u_norm = VecNorm(u);
    if( fabs(u_norm-1.0) > 1e-5 ){
        if( p_rank == 0 ){
            printf("Error in GenRotationMatrix: u_norm differs too much from 1.0!\n");
        }
    }
    
    R[0] = cos_theta+ux*ux*One_cos_theta;
    R[1] = ux*uy*One_cos_theta - uz*sin_theta;
    R[2] = ux*uz*One_cos_theta + uy*sin_theta;

    R[3] = uy*ux*One_cos_theta + uz*sin_theta;
    R[4] = cos_theta + uy*uy*One_cos_theta;
    R[5] = uy*uz*One_cos_theta - ux*sin_theta;

    R[6] = uz*ux*One_cos_theta - uy*sin_theta;
    R[7] = uz*uy*One_cos_theta + ux*sin_theta;
    R[8] = cos_theta + uz*uz*One_cos_theta;

    return 0;
}

//  这个版本在“新旧”指向的经度差值大于90度的时候给出的转动角不是最优化的结果。
//  该版本已经被弃用！
int Get_RotationAngle_faster( double ra_old, double dec_old,
                              double ra_new, double dec_new,
                              double* angle_rot ){
//  所有角度都是以“弧度”为单位的！！！
    int info = 0;
    double sin_ra_old = sin(ra_old*PI_180);
    double cos_ra_old = cos(ra_old*PI_180);
    double sin_dec_old = sin(dec_old*PI_180);
    double cos_dec_old = cos(dec_old*PI_180);

    double sin_ra_new = sin(ra_new*PI_180);
    double cos_ra_new = cos(ra_new*PI_180);
    double sin_dec_new = sin(dec_new*PI_180);
    double cos_dec_new = cos(dec_new*PI_180);

//  compute the rotation angle
    double Trace_of_R = cos_dec_old*cos_ra_old*cos_ra_new*cos_dec_new
                      + sin_dec_old*sin_dec_new
                      + cos_dec_old*sin_ra_old*sin_ra_new*cos_dec_new
                      + sin_dec_old*cos_ra_old*cos_ra_new*sin_dec_new
                      + cos_dec_new*cos_dec_old
                      + sin_dec_old*sin_ra_old*sin_ra_new*sin_dec_new
                      + sin_ra_old*sin_ra_new
                      + cos_ra_old*cos_ra_new;

    double cosval     = 0.5*(Trace_of_R-1.0);

    cosval = (cosval >  1.0) ?  1.0 : cosval;
    cosval = (cosval < -1.0) ? -1.0 : cosval;

    *angle_rot = acos(cosval)/PI_180;

    return info;
}

//  ==================================================================
//  原先所采取的旋转方式在某些情况下给出的旋转角度不是最合理的。例如在靠近黄极的时候，
//  如果两个指向的经度相差接近于180度，利用原先的旋转方式得到旋转角度就接近于180度，
//  但实际上可以通过更小角度的转动来完成指向的移动。
//
//  注意：新版本返回的角度单位是“度数”,不再使用“弧度”！
//
//  另外，目前还没有优化这个函数。
int Get_RotationAngle_faster2( double ra_old, double dec_old,
                               double ra_new, double dec_new,
                               double* angle_rot,
                               int rank ){

    if( gsl_isnan(ra_old) || gsl_isnan(dec_old) ){
        if( rank == 0 ){
            printf("$$$$$$$$$$\n");
        }
        MPI_Finalize();
        exit(0);
    }

    int info = 0;

    double sin_ra_old = sin(ra_old*PI_180);
    double cos_ra_old = cos(ra_old*PI_180);
    double sin_dec_old = sin(dec_old*PI_180);
    double cos_dec_old = cos(dec_old*PI_180);

    double sin_ra_new = sin(ra_new*PI_180);
    double cos_ra_new = cos(ra_new*PI_180);
    double sin_dec_new = sin(dec_new*PI_180);
    double cos_dec_new = cos(dec_new*PI_180);

    double p_old[3] = {sin_dec_old*cos_ra_old,sin_dec_old*sin_ra_old,cos_dec_old};
    double p_new[3] = {sin_dec_new*cos_ra_new,sin_dec_new*sin_ra_new,cos_dec_new};
    double p_tmp[3] = {-999,-999,-999};

// #if defined(_USE_OLD_ROTATION_)
//     double cosval = p_old[0]*p_new[0] + p_old[1]*p_new[1] + p_old[2]*p_new[2];
//     if(cosval>1.0) cosval = 1.0;
//     if(cosval<-1.0) cosval = -1.0;

//     *angle_rot = acos(cosval)*180/PI;
//     return info;
// #endif

    double Rz[9] = {0};  //  给存贮绕z-轴的两个旋转矩阵分配内存

//  ====================================================================
//  step 1:绕z轴旋转，将旧指向旋转到新指向所在的子午圈，并计算出旋转后的“临时”指向矢量

    double delta_alpha=-999;
    double delta_alpha1,delta_alpha2;

    if( fabs(ra_new-ra_old) < 1e-5 ){ //  这个情况下可以近似认为两个指向处在同一个子午圈内

        if( fabs(dec_new-dec_old) < 1e-5 ){
            *angle_rot = 0;
            return info; //两个指向相同的情况下不需要进行任何转动。
        }

        delta_alpha = 0;
        delta_alpha1 = 0;   //等价于不做任何转动。
        delta_alpha2 = 180; //这种情况应该予以排除，因为会导致帆板面的指向反转，造成无法接受太阳光照进行发电。

        p_tmp[0] = p_old[0];
        p_tmp[1] = p_old[1];
        p_tmp[2] = p_old[2];

        double alpha_n[3] = {0,0,1};
        GenRotationMatrix(alpha_n,0,Rz,rank);
    }
    else{
        //  alpha_old,alpha_new是在XY平面内的单位向量
        double alpha_old[3] = {cos_ra_old, sin_ra_old, 0};
        double alpha_new[3] = {cos_ra_new, sin_ra_new, 0};
        
        //  根据alpha_old,alpha_new来确定如何绕z-轴旋转
        double cosval = DotProduct(alpha_old,alpha_new);
        cosval = (cosval >  1.0) ?  1.0 : cosval;
        cosval = (cosval < -1.0) ? -1.0 : cosval;
        delta_alpha1 = acos(cosval)*180/PI; //这个值始终是大于0的
        delta_alpha2 = 180-delta_alpha1; //这个值始终也是大于0的

        double alpha_n[3] = {-999,-999,-999};
        CrossProduct(alpha_old,alpha_new,alpha_n);

        double alpha_n_norm = VecNorm(alpha_n);
        alpha_n[0] /= alpha_n_norm;
        alpha_n[1] /= alpha_n_norm;
        alpha_n[2] /= alpha_n_norm;

    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  因为旋转角超过90度时帆板面的指向会“反转”，导致不能接受太阳光照，因此只考虑
    //  旋转角小于90度的情况。
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if( fabs(delta_alpha1) <= 90 ){
            delta_alpha = delta_alpha1;
            // printf("---> using delta_alpha1 = %g (deg) to generate the rotation matrix ...\n",delta_alpha1);
            // print_vec("alpha_n",alpha_n);
            //转动的方向由转动轴的指向决定，因此传入的角度值都是正数
            GenRotationMatrix(alpha_n,delta_alpha1,Rz,rank);
        }
        else if( fabs(delta_alpha2) <= 90 ){
            delta_alpha = delta_alpha2;
            // printf("---> using delta_alpha2 = %g (deg) to generate the rotation matrix ...\n",delta_alpha2);
            //转动的方向由转动轴的指向决定，因此传入的角度值都是正数
            //这里需要注意的是需要将alpha_n反号
            alpha_n[0] *= -1;
            alpha_n[1] *= -1;
            alpha_n[2] *= -1;
            GenRotationMatrix(alpha_n,delta_alpha2,Rz,rank);
        }
        else{
            ErrorInfo("failed to calculation rotation matrix Rz");
            printf("ra_old = %.10g, dec_old = %.10g\n",ra_old,dec_old);
            printf("ra_new = %.10g, dec_new = %.10g\n",ra_new,dec_new);
            printf("delta_alpha1 = %g\n",delta_alpha1);
            printf("delta_alpha2 = %g\n",delta_alpha2);
            MPI_Finalize();
            exit(0);
        }
    
        MatrixVecProduct(Rz,p_old,p_tmp);
    }

//  ====================================================================
//  step 2:根据alpha_tmp和alpha_new的叉乘来计算第二次旋转的旋转轴；虽然确定了旋转轴，
//  但是依旧分别存在两种旋转方式。

    double axis[3] = {0,0,0};

    //  首先还是先判断“临时指向”与目标指向是否共线(包括同向与反向)
    double p_tmp_dot_p_new = DotProduct(p_tmp,p_new);

    if( fabs(p_tmp_dot_p_new-1) < 1e-10 ){// p_tmp与p_new重合
        *angle_rot = delta_alpha;
        return info;
    }
    else if( fabs(p_tmp_dot_p_new+1) < 1e-10 ){// p_tmp与p_new反向
        //  旋转轴位于XY平面内，因此可以直接根据ra_new+pi/2直接计算出来旋转轴；此时
        //  旋转轴的指向不影响结果
        axis[0] = cos(0.5*PI+ra_new*PI_180);
        axis[1] = sin(0.5*PI+ra_new*PI_180);
        axis[2] = 0;

        double R[9] = {0};
        double Rn[9] = {0};

        GenRotationMatrix(axis,180,Rn,rank);
        MatrixMultiplication(Rn,Rz,R);  //  此处需要注意旋转矩阵乘积的顺序
        *angle_rot = GetRotationAngleFromMatrix(R);

        return info;
    }
    else{
        //  给存贮绕n-轴的两个旋转矩阵分配内存（n轴由两个矢量的叉积得到）
        double R1[9] = {0};
        // double R2[9] = {0};
        double RAll1[9] = {0};
        // double RAll2[9] = {0};

        double cosval = DotProduct(p_tmp,p_new);
        cosval = (cosval >  1.0) ?  1.0 : cosval;
        cosval = (cosval < -1.0) ? -1.0 : cosval;
        delta_alpha1 = acos(cosval)*180/PI; //这个值始终是大于0的
        delta_alpha2 = 360-delta_alpha1; //这个值始终也是大于0的，但表示反方向旋转

        CrossProduct(p_tmp,p_new,axis);
        double axis_norm = VecNorm(axis);
        axis[0] /= axis_norm;
        axis[1] /= axis_norm;
        axis[2] /= axis_norm;

        GenRotationMatrix(axis,delta_alpha1,R1,rank);
        MatrixMultiplication(R1,Rz,RAll1);

        *angle_rot = GetRotationAngleFromMatrix(RAll1);
        return info;
    }

    return info;
}

int Get_RotationAngle_Zhang( double ra_old, double dec_old, double ra_new, double dec_new, double *angle ){

	double ppNorm[3],p1[3],p2[3],pn1[3],pn2[3];
	p1[0] = sin(dec_old*PI_180)*cos(ra_old*PI_180);
	p1[1] = sin(dec_old*PI_180)*sin(ra_old*PI_180);
	p1[2] = cos(dec_old*PI_180);

	p2[0] = sin(dec_new*PI_180)*cos(ra_new*PI_180);
	p2[1] = sin(dec_new*PI_180)*sin(ra_new*PI_180);
	p2[2] = cos(dec_new*PI_180);

	crossMultiple(p1,p2,ppNorm);
	crossMultiple(p1,ppNorm,pn1);

	crossMultiple(p2,ppNorm,pn2);
	*angle = getAngle132(pn1[0],pn1[1],pn1[2],pn2[0],pn2[1],pn2[2],0,0,0);

	return 0;
}


//  按照给定的新旧指向，生成一个将旧指向往新指向转动5度的旋转矩阵
int GenRotationMatrix5deg(  double ra_old, double dec_old, 
                            double ra_new, double dec_new,
                            double R5deg[],
                            double rot_axis[],
                            int p_rank ){

    double sin_ra_old = sin(ra_old*PI_180);
    double cos_ra_old = cos(ra_old*PI_180);
    double sin_dec_old = sin(dec_old*PI_180);
    double cos_dec_old = cos(dec_old*PI_180);

    double sin_ra_new = sin(ra_new*PI_180);
    double cos_ra_new = cos(ra_new*PI_180);
    double sin_dec_new = sin(dec_new*PI_180);
    double cos_dec_new = cos(dec_new*PI_180);

    double p_old[3] = { cos_dec_old*cos_ra_old,
                        cos_dec_old*sin_ra_old,
                        sin_dec_old };
    double p_new[3] = { cos_dec_new*cos_ra_new,
                        cos_dec_new*sin_ra_new,
                        sin_dec_new };
    double p_tmp[3] = {-999,-999,-999};

    double Rz[9] = {0};  //  给存贮绕z-轴的两个旋转矩阵分配内存
    double RAll[9] = {0};
    double cosval;

//  ====================================================================
//  step 1:绕z轴旋转，将旧指向旋转到新指向所在的子午圈，并计算出旋转后的“临时”指向矢量

    double delta_alpha1,delta_alpha2;

    if( fabs(ra_new-ra_old) < 1e-5 ){ //  这个情况下可以近似认为两个指向处在同一个子午圈内

        if( fabs(dec_new-dec_old) < 1e-5 ){  //两个指向相同的情况下不需要进行任何转动
            RAll[0] = 1;
            RAll[1] = 0;
            RAll[2] = 0;

            RAll[3] = 0;
            RAll[4] = 1;
            RAll[5] = 0;
            
            RAll[6] = 0;
            RAll[7] = 0;
            RAll[8] = 1;

            R5deg[0] = 1;
            R5deg[1] = 0;
            R5deg[2] = 0;
            
            R5deg[3] = 0;
            R5deg[4] = 1;
            R5deg[5] = 0;
            
            R5deg[6] = 0;
            R5deg[7] = 0;
            R5deg[8] = 1;

            return 0;
        }

        p_tmp[0] = p_old[0];
        p_tmp[1] = p_old[1];
        p_tmp[2] = p_old[2];
        double alpha_n[3] = {0,0,1};
        GenRotationMatrix(alpha_n,0,Rz,p_rank); // Rz is actually a 3x3 unit matrix
    }
    else{
        //  alpha_old,alpha_new是在XY平面内的单位向量
        double alpha_old[3] = {cos_ra_old, sin_ra_old, 0};
        double alpha_new[3] = {cos_ra_new, sin_ra_new, 0};
        
        //  根据alpha_old,alpha_new来确定如何绕z-轴旋转
        cosval = DotProduct(alpha_old,alpha_new);
        cosval = (cosval >  1.0) ?  1.0 : cosval;
        cosval = (cosval < -1.0) ? -1.0 : cosval;
        delta_alpha1 = acos(cosval)*180/PI; //这个值始终是大于0的
        delta_alpha2 = 180-delta_alpha1; //这个值始终也是大于0的

        double alpha_n[3] = {-999,-999,-999};
        CrossProduct(alpha_old,alpha_new,alpha_n);

        double alpha_n_norm = VecNorm(alpha_n);
        alpha_n[0] /= alpha_n_norm;
        alpha_n[1] /= alpha_n_norm;
        alpha_n[2] /= alpha_n_norm;

    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  因为旋转角超过90度时帆板面的指向会“反转”，导致不能接受太阳光照，因此只考虑
    //  旋转角小于90度的情况。
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if( fabs(delta_alpha1) <= 90 ){
            GenRotationMatrix(alpha_n,delta_alpha1,Rz,p_rank);
        }
        else if( fabs(delta_alpha2) <= 90 ){
            alpha_n[0] *= -1;
            alpha_n[1] *= -1;
            alpha_n[2] *= -1;
            GenRotationMatrix(alpha_n,delta_alpha2,Rz,p_rank);
        }
        else{
        //  也许由于数值原因，导致失败，此时返回失败码（非0即可，后期debug的时候可以根据错误码进行定位）
            print_Error_Msg("error code -111, failed to compute rotation matrix Rz[]",p_rank);
            return -111;
        }
    
        MatrixVecProduct(Rz,p_old,p_tmp);
    }

//  ====================================================================
//  step 2:根据alpha_tmp和alpha_new的叉乘来计算第二次旋转的旋转轴；虽然确定了旋转轴，
//  但是依旧分别存在两种旋转方式。

    double axis[3] = {0,0,0};
    double Rn[9] = {0};

    // print_vec("p_old",p_old);
    // print_vec("p_new",p_new);
    double p_tmp_dot_p_new = DotProduct(p_tmp,p_new);
    // printf("p_tmp_dot_p_new = %20.15f\n",p_tmp_dot_p_new);

    //  首先判断“临时指向”与目标指向是否共线(包括同向与反向)
    if( fabs(p_tmp_dot_p_new-1) < 1e-10 ){// p_tmp与p_new重合
        RAll[0] = Rz[0];
        RAll[1] = Rz[1];
        RAll[2] = Rz[2];
        RAll[3] = Rz[3];
        RAll[4] = Rz[4];
        RAll[5] = Rz[5];
        RAll[6] = Rz[6];
        RAll[7] = Rz[7];
        RAll[8] = Rz[8];
    }
    else if( fabs(p_tmp_dot_p_new+1) < 1e-10 ){// p_tmp与p_new反向
        //  旋转轴位于XY平面内，因此可以直接根据ra_new+pi/2直接计算出来旋转轴；此时
        //  旋转轴的指向不影响结果
        axis[0] = cos(0.5*PI+ra_new*PI_180);
        axis[1] = sin(0.5*PI+ra_new*PI_180);
        axis[2] = 0;

        GenRotationMatrix(axis,180,Rn,p_rank);
        MatrixMultiplication(Rn,Rz,RAll);  //  此处需要注意旋转矩阵乘积的顺序
    }
    else{
        cosval = DotProduct(p_tmp,p_new);
        cosval = (cosval >  1.0) ?  1.0 : cosval;
        cosval = (cosval < -1.0) ? -1.0 : cosval;
        delta_alpha1 = acos(cosval)*180/PI; //这个值始终是大于0的
        delta_alpha2 = 360-delta_alpha1; //这个值始终也是大于0的，但表示反方向旋转

        if( gsl_isnan(delta_alpha1) || gsl_isnan(delta_alpha2) ){
            print_Error_Msg("error code -222, either delta_alpha1 or delta_alpha2 is NAN",p_rank);
            return -222;
        }

        CrossProduct(p_tmp,p_new,axis);
        double axis_norm = VecNorm(axis);
        axis[0] /= axis_norm;
        axis[1] /= axis_norm;
        axis[2] /= axis_norm;

        GenRotationMatrix(axis,delta_alpha1,Rn,p_rank);
        MatrixMultiplication(Rn,Rz,RAll);
    }

    double u[3] = {0,0,0};
    get_rotation_axis(RAll,u);  //  获取旋转轴
    double u_norm = VecNorm(u); //  归一化旋转轴

    if( gsl_isnan(u_norm) == 1 || fabs(u_norm) < 1e-10){
        print_Error_Msg("error code -333, failed to get the rotation axis u[3]",p_rank);
        printf("%.20g\n",u_norm);
        return -333;
    }

    u[0] /= u_norm;
    u[1] /= u_norm;
    u[2] /= u_norm;

    rot_axis[0] = u[0];
    rot_axis[1] = u[1];
    rot_axis[2] = u[2];

    GenRotationMatrix( u, 5.0, R5deg, p_rank );  //  生成相应的旋转矩阵：
    return 0;
}

//  将望远镜的指向按照旋转矩阵R5deg[] 旋转5度
int Rotate_about_AxisP_by_5deg( double R5deg[], 
                                double ra_old, 
                                double dec_old, 
                                double *ra_new, 
                                double *dec_new,
                                int p_rank ){
    double radec[2];

    double sin_ra_old = sin(ra_old*PI_180);
    double cos_ra_old = cos(ra_old*PI_180);
    double sin_dec_old = sin(dec_old*PI_180);
    double cos_dec_old = cos(dec_old*PI_180);

    double x_old[3] = { cos_dec_old*cos_ra_old,
                        cos_dec_old*sin_ra_old,
                        sin_dec_old };

    double x_new[3] = {0,0,0};

    MatrixVecProduct(R5deg,x_old,x_new);
    Cartesian2Equatorial(x_new,radec);
    *ra_new = radec[0];
    *dec_new= radec[1];

    if( gsl_isnan(radec[0]) || gsl_isnan(radec[1]) ){
        if( p_rank == 0 ){
            print_Error_Msg(" either radec[0] or radec[1] is NAN",p_rank);
            print_mat("==> R5deg", R5deg);
            print_vec("==> x_old", x_old);
            print_vec("==> x_new", x_new);
        }
    }

    // print_mat("==> R5deg", R5deg);
    // print_vec("==> x_old", x_old);
    // print_vec("==> x_new", x_new);

    // printf("# ra_new = %.10g, dec_new = %.10g\n",*ra_new,*dec_new);

    return 0;
}
