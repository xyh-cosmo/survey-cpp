#!/bin/bash

#qsub -q i7.q -pe mpich 12 0720_G17_E15.job  
#qsub -q i8.q -pe mpich 12 0720_G17_E16.job  
#qsub -q i9.q -pe mpich 12 0720_G17_E17.job
#qsub -q i0.q -pe mpich 12 0720_G17_E18.job  
#qsub -q i1.q -pe mpich 12 0720_G17_E19.job  
#qsub -q i2.q -pe mpich 12 0720_G17_E20.job


qsub -q i1.q -pe mpich 12 0720_G17_E20.jobX_0yr
qsub -q i2.q -pe mpich 12 0720_G17_E20.jobX_4yr
qsub -q i3.q -pe mpich 12 0720_G17_E20.jobX_5yr
qsub -q i4.q -pe mpich 12 0720_G17_E20.jobX_6yr
qsub -q i5.q -pe mpich 12 0720_G17_E20.jobX_7yr
qsub -q i6.q -pe mpich 12 0720_G17_E20.jobX_8yr


#qsub -q i7.q -pe mpich 12 0720_G17_E17.jobY
#qsub -q i8.q -pe mpich 12 0720_G17_E18.jobY
#qsub -q i9.q -pe mpich 12 0720_G17_E19.jobY
#qsub -q i0.q -pe mpich 12 0720_G17_E20.jobY
