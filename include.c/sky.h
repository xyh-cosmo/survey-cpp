#ifndef _SKY_H_
#define _SKY_H_

typedef enum {
	deep,
	ultraDeep,
	deepAndUltra,
	spectal
} SurveyType;

struct Coordinate_SM {

//	最后选择权重最小的作为最优的选择结果（不是权重最大的）
	double weight;  //  计算权重，这两个可以合起来
	double wFactor; //  计算权重，这两个可以合起来

	double dec; // （黄）纬度
	double ra;	// （黄）经度
	double inDeepFlag; //是否在极深度成像区域
	double elat; // ecliptic latitude， 黄纬度，之前用赤经赤纬坐标，需要这个数，目前这个数和dec是一个数值
	double gb; // galaxy latitude， 银纬

	//Cartesian coordinate ,related with equator coordinate，笛卡尔坐标，第四个是距离原点（也就是地心）的距离
	double x;
	double y;
	double z;
	double dist2orig;
	/*-------------------solar panel-------------------*/
	//	normal vector of the plane of solar panel
	// 	这三个是帆板法线
	double nx;
	double ny;
	double nz;
	//	normal vector +-30 deg of the plane of solar panel
	// 	这六个分别是帆板法线在法线平面内向两个方向转动一定角度后的向量，之前是转30°，现在不一定是30°
	double nx30;
	double ny30;
	double nz30;
	double nx30_;
	double ny30_;
	double nz30_;
	//the normal of the plane(decided by the vector point and the normal of solar panel)
	// 有这样一个平面，指向向量和帆板法线向量构成的平面，下面三个是这面的法线
	double nnx;
	double nny;
	double nnz;
	/*-------------------END solar panel-------------------*/

	//两个方向的间隔，可能不用了
	double interval_lon;
	double interval_lat;

	int flag;	            //标志位，目前用于记录该天区观测次数
	int id;	                //每一个天区都有一个id，从0开始编号的，前面多色成像天区，后面极深度
	int IsInSunSideFlag;    // 是否在阳照区
    int isObserve;	        // 0表示不可观测； 1表示可观测
	int targetCoverNum;	    //该天区需要观测次数
	int maxCoverNum;	    //该天区最大可以观测次数

	// 下面四个是该天区上下左右相邻的四个天区的id，在高纬度不是很准确，高纬度找最近的四个天区
	int left_neighbour;
	int right_neighbour;
	int up_neighbour;
	int down_neighbour;

            //0   1  2  3  4  5  6  7    8    9
//10 filter, NUV, u ,g ,r ,i ,z ,y ,GI , GU , GV
	int filter[10]; //记录被滤光片观测次数,应该是只是才策略上用到了，真正计算每个区域被多少滤光片覆盖用的不是这个，用的是HealPix那部分的，这个在高纬不准确

} CoordSM;

typedef struct Coordinate_SM SKY_Coord;

void init_SKY_Coord( SKY_Coord* sky );

int* readData(char*, int, int*);
void copySkyMapData(SKY_Coord*, SKY_Coord*);
SKY_Coord* produceSkyArea(int*);
int* readData(char* fileName, int elemNum, int* datNum);
void copySkyMapData(SKY_Coord* from, SKY_Coord* to);
void calculateSkySize(double size[2],double x_size, double y_size,double x_overLap, double y_overlap);
int SplitSkyArea( int yLen, 
                  int xLen, 
                  SKY_Coord* skyMap,
                  double x_size, 
                  double y_size );
SKY_Coord* produceSkyArea(int* normal_sky_num);

//	每隔一段时间更新可观测天区的id
int sky_num_remained;	//对剩下的可观测天区的数目进行计数
int* sky_id_tracker;	// 用于保存可观测天区的ID
int* update_sky_id_tracker( SKY_Coord* sky, int sky_num, int rank, int* remained_num );

int count_sky_num( SKY_Coord* sky, int sky_num, int rank );

#endif