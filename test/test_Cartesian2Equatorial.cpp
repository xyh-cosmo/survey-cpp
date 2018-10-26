//  Purpose of this test:
//  The original transformation seems uneffcient, therefore I suspect it can be replaced by 
//  a new and much faster one!
//  
//  Test Date: Mar-6,2018
//
//  Test details:
//  Using randomly generated 10000 Cartesian coordiantes, and then transform these coordiantes
//  into Equatorial coordinates (3D->2D), and compare the time costs as well as the precisions.

#include <iostream>
#include <cmath>
#include <random>

using namespace std;

#define PI 3.1415926535
#define N_TOT 10000

void Cartesian2Equatorial(double* carCoor, double* eCoor) {
    if (carCoor[0] > 0 && carCoor[1] >= 0) {
        *eCoor = atan(carCoor[1] / carCoor[0]) * 360 / (2 * PI);
    } else if (carCoor[0] < 0) {
        *eCoor = (atan(carCoor[1] / carCoor[0]) + PI) * 360 / (2 * PI);
    } else if (carCoor[0] > 0 && carCoor[1] < 0) {
        *eCoor = (atan(carCoor[1] / carCoor[0]) + 2 * PI) * 360 / (2 * PI);
    } else if (carCoor[0] == 0 && carCoor[1] < 0) {
        *eCoor = 270;
    } else if (carCoor[0] == 0 && carCoor[1] > 0) {
        *eCoor = 90;
    }
    *(eCoor + 1) = atan(carCoor[2] / sqrt(carCoor[0] * carCoor[0] + carCoor[1] * carCoor[1])) * 360 / (2 * PI);
}


void NewCar2Eql(double* carCoor, double* eCoor){
    double x1 = carCoor[0];
    double x2 = carCoor[1];
    double x3 = carCoor[2];
    
    double r = sqrt(x1*x1+x2*x2+x3*x3);
    
    double theta = asin(x3/r);
    *(eCoor+1) = theta*180./PI;
    *(eCoor+0) = atan(x2/(r*cos(theta)+x1)) *360./PI;
    
    if( *(eCoor+0) < 0 )
        *(eCoor+0) += 360.;
}

int main(int argc, char *argv[]){

    unsigned int seed = time(NULL);
    default_random_engine dre(seed);
    uniform_real_distribution<double> dr(-100.,100.);

    
    double t_start, t_end;
    double time_cost1, time_cost2;
    
    cout << "preparing random Cartesian coordinates ...\n";
    double* x1 = new double[N_TOT];
    double* x2 = new double[N_TOT];
    double* x3 = new double[N_TOT];
    
    for( int i=0; i<N_TOT; ++i ){
        x1[i] = dr(dre);
        x2[i] = dr(dre);
        x3[i] = dr(dre);
    }
    
    cout << "done!\n";

    double Ra1[N_TOT], Dec1[N_TOT];
    double Ra2[N_TOT], Dec2[N_TOT];
    double carCoor[3], eCoor[2];
    
    t_start = time(NULL);
    for( int i=0; i<N_TOT; ++i ){
        carCoor[0] = x1[i];
        carCoor[1] = x2[i];
        carCoor[2] = x3[i];

//  =================================================        
        Cartesian2Equatorial(carCoor,eCoor);
        Ra1[i] = eCoor[0];
        Dec1[i] = eCoor[1];
    }
    t_end = time(NULL);
    time_cost1 = t_end - t_start;

    t_start = time(NULL);
    for( int i=0; i<N_TOT; ++i ){
        carCoor[0] = x1[i];
        carCoor[1] = x2[i];
        carCoor[2] = x3[i];
        
        NewCar2Eql(carCoor,eCoor);
        Ra2[i] = eCoor[0];
        Dec2[i] = eCoor[1];
    }
    t_end = time(NULL);
    time_cost2 = t_end - t_start;


//    for( int i=0; i<N_TOT; ++i ){
//        cout << "err_Ra = " << Ra2[i]-Ra1[i]
//             << "\terr_Dec = " << Dec2[i]-Dec1[i] << endl;
//    }

    cout << "time_cost1 = " << time_cost1 << endl;
    cout << "time_cost2 = " << time_cost2 << endl;
    
    delete[] x1;
    delete[] x2;
    delete[] x3;
    
    return 0;
}




