#   ============================================================================
#   Youhua Xu Dec-13-2017
#   ============================================================================
#   Make a build dir for compilation
MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
BINDIR = $(MDIR)/bin
DEBUG  = $(MDIR)/debug

#	make build & binary dirs
.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; fi;
	touch $(WRKDIR)/.base;
	if ! [ -e $(BINDIR) ]; then mkdir $(BINDIR) ; fi;
	touch $(BINDIR)/.base;
	if ! [ -e $(DEBUG) ]; then mkdir $(DEBUG) ; fi;
	touch $(DEBUG)/.base;
#   ============================================================================
#   Set the source file path
vpath %.cpp main:src:test
vpath %.cc  main:src:test
vpath %.c   main:src:test
vpath %.hpp include
vpath %.h   include
vpath %.o build
vpath .base build

INCLUDES = -I $(MDIR)/include

CPP			= gcc
CC          = mpicc
CCFLAG  	= -Wall -DHAVE_INLINE -D_FULL_DEBUG_
CCFLAG += -D_USE_SWITCH_
# CCFLAG += -pg
CCFLAG	+= -D_ENABLE_BETA_ANGLE_  # 控制是否开启beta角限制

#CCFLAG	+= -D_USE_OLD_SKY_SPLIT_  # 控制天区划分的方式
CCFLAG	+= -D_USE_SKY_SPLIT_2_

CCFLAG 	+= -D_ENABLE_HIGH_LATITUDE_PRIOR_ # 控制是否开启高纬度优先的策略
#CCFLAG	+= -D_ENABLE_MARK_OBS_POINT_  # 控制是否开启实时的面积统计
# CCFLAG 	+= -D_ENABLE_GRI_WEIGHT_  # 根据gri波段的覆盖情况来修改权重
# CCFLAG	+= -D_ENABLE_DEEP_PRIOR_  # 控制是否编译“极深场优先”的那部分代码
CCFLAG += -D_ENABLE_CONTINOUS_OBS_  # 控制是否开启“优先观测相邻天区”的策略
CCFLAG += -D_OLD_CMG_USE_  # 控制是否使用张鑫的旧版本CMG use的计算
#CCFLAG += -D_USE_OLD_ROTATION_
CCFLAG += -D_UPDATE_ONLY_SKY_NUM_
CCFLAG += -D_ENABLE_SAVE_SKYMAP_
CCFLAG += -D_ENABLE_PANEL_ROTATION_ # 控制是否开启帆板转动
CCFLAG += -D_TURN_ON_5DEG_DRIFT_ #控制是否开启5度转动“牵引”，用于测试FindNewPointByAllSearch3
CCFLAG += -D_TURN_DEC60_PRIOR_
#CCFLAG += -D_ENABLE_DARKSIDE_PANEL_SUN_ANGLE_CHECK_
#CCFLAG += -D_USE_CMG_ONE_ORBIT_STATE_
#CCFLAG += -D_SLEW_TIME_GSL_INTERP_  # 控制是否使用新版本的转动时间函数

# CCFLAG	+= -D_USE_LAGRANGE_INTERP_

### 天区搜索函数的版本控制：
# CCFLAG += -D_ENABLE_NEW_SEARCH_	# 控制是否使用天区搜索函数（FindNewPointByAllSearch）
# CCFLAG += -D_ENABLE_NEW_SEARCH2_	# 控制是否使用新的天区搜索函数（FindNewPointByAllSearch2）
CCFLAG += -D_ENABLE_NEW_SEARCH3_	# 控制是否使用新的天区搜索函数（FindNewPointByAllSearch3）
#CCFLAG += -D_ENABLE_NEW_SEARCH4_	# 控制是否使用新的天区搜索函数（FindNewPointByAllSearch4）

OPTFLAG		= -O2 #-ffast-math #( not recognized by intel compiler )
# OPTFLAG += -pg

LDFLAG      =
#   http://www.tuicool.com/articles/NBfeYj
ifeq ($(shell uname -s),Linux)
	LDFLAG	+= -Wl,--no-as-needed
endif
#LIBS		= -limcmc
#LIBS		+= -larmadillo -lblas -llapack
#LIBS        += -lgsl -lgslcblas

LIBS	= -lgsl -lgslcblas -lm

ifeq ($(shell uname -s),Darwin)
	LDFLAG	+= -framework Accelerate #(-framework Accelerate is for Mac OSX)
endif

%.o: %.cpp .base
	cd $(WRKDIR);$(CC) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.o: %.cc .base
	cd $(WRKDIR);$(CC) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.o: %.c .base
	cd $(WRKDIR);$(CC) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

VERBOSE = Verbose.o

SURVEY_SIM_MPI	= SurveySim_MPI.o
COVERSTAT       = CoverStat.o
COVERSTAT_A     = CoverStat_A.o
COVERSTAT_B     = CoverStat_B.o
COVERSTAT_C     = CoverStat_C.o
GET_BETA_TIME   = get_beta_time.o
ADJ_BETA_TIME	= adjust_beta.o
EXTRACT_DT			= ExtractTimeInterval.o

# useful tools
CALC_ANGLE = CalcAngle.o
GALACTIC2ECLIPTIC= Galactic2Ecliptic.o
ECLIPTIC2GALACTIC= Ecliptic2Galactic.o

#  tests and debugs
CHECK_ROTATION	= check_rotation.o
CMG_CHECK		= CMG_check.o
TEST_ANGLE_CALC = TestAngleCalc.o
RAND_TEST_ANGLE_CALC = RandTestAngleCalc.o
DEBUG_PANEL_ANGLE = DebugTestPanelAngle.o
GET_TRANS_TIME = GetTransTime.o
TEST_GENROTATIONMATRIX5DEG = test_GenRotationMatrix5deg.o
DEBUG_JPL405 = TestJPL405.o
TEST_ORBIT_TIME_COST = TestSatOrbitTimeCost.o
TEST_ORBIT_TIME_COST2 = TestSatOrbitTimeCost2.o
CHECK_CART2SPH = Cart2Spherical.o
DEBUG_OBS_MOON = DebugIsObscureByMoon.o
TEST_UNDER_SAT = TestUnderSat.o
TEST_DRIFT = TestDrift.o

DEPS = LocateSatellite.o transformTools.o SkyAreaSplit.o SurveyConditionLimit.o \
		init.o extra.o Parser.o CMG.o CoverAreaStat.o CoverAreaStat2.o  rotation.o \
		ccd.o TestConditions.o Panel.o drift5deg.o satellite.o \
		SearchNewPoints.o SearchNewPoints_new2.o SearchNewPoints_new3.o SearchNewPoints_new4.o
EPHCMP = gnulliver.o ephcom.o

# all:SurveySim GetTime adjust_beta check_rotation cover_stat_a cover_stat_b cover_stat_c \
# 	test_angle_calc rand_test_angle_calc extract_dT calcAngle cmg_check debug_panel_angle \
# 	get_tTime g2el el2g test_GenRotationMatrix5deg #cover_stat

all: tools debug

tools:SurveySim get_beta_time adjust_beta cover_stat_a cover_stat_b cover_stat_c \
	extract_dT calcAngle get_tTime g2el el2g

debug:check_rotation test_angle_calc rand_test_angle_calc cmg_check debug_panel_angle \
	test_GenRotationMatrix5deg test_orbit_time_cost test_orbit_time_cost2 cart2shp \
	debug_osb_moon test_under_sat test_drift #cover_stat


SurveySim:${SURVEY_SIM_MPI} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

extract_dT:${EXTRACT_DT} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

cover_stat:${COVERSTAT} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

cover_stat_a:${COVERSTAT_A} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

cover_stat_b:${COVERSTAT_B} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

cover_stat_c:${COVERSTAT_C} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

get_beta_time:${GET_BETA_TIME} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

get_tTime:${GET_TRANS_TIME} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

adjust_beta:${ADJ_BETA_TIME} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

g2el:${GALACTIC2ECLIPTIC} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@

el2g:${ECLIPTIC2GALACTIC} ${EPHCMP} ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@
#   ================================================================================================
# 	The following are for debugging!
cart2shp:${CHECK_CART2SPH}  ${EPHCMP}  ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

check_rotation:${CHECK_ROTATION}  ${EPHCMP}  ${DEPS}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

test_angle_calc:${TEST_ANGLE_CALC} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

rand_test_angle_calc:${RAND_TEST_ANGLE_CALC} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

calcAngle:${CALC_ANGLE} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@
	
cmg_check:${CMG_CHECK} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

debug_panel_angle:${DEBUG_PANEL_ANGLE} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

test_GenRotationMatrix5deg:${TEST_GENROTATIONMATRIX5DEG} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

debug_jpl405:${DEBUG_JPL405} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

test_orbit_time_cost:${TEST_ORBIT_TIME_COST} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

test_orbit_time_cost2:${TEST_ORBIT_TIME_COST2} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

test_under_sat:${TEST_UNDER_SAT} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

debug_osb_moon:${DEBUG_OBS_MOON} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@

test_drift:${TEST_DRIFT} ${DEPS} ${EPHCMP}
	${CC} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(DEBUG)/$@
#   ================================================================================================
.PHONY:clean tidy run
clean: .base
	rm -rf $(WRKDIR);
tidy:
	make clean; rm -rf $(BINDIR); rm -rf $(DEBUG)
