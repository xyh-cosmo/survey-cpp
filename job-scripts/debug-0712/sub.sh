#!/bin/bash

qsub -q i0.q -pe mpich 12 B18_E19.5.job
#qsub -q i1.q -pe mpich 12 B18_E20.job
qsub -q i1.q -pe mpich 12 B19_E19.5.job

qsub -q i2.q -pe mpich 12 B18_E20.5.job
qsub -q i3.q -pe mpich 12 B18_E21.job

#qsub -q i4.q -pe mpich 12 B19_E18.5.job
#qsub -q i5.q -pe mpich 12 B19_E19.job
#qsub -q i6.q -pe mpich 12 B19_E19.5.job
#qsub -q i7.q -pe mpich 12 B19_E20.job

qsub -q i4.q -pe mpich 12 B18_E19.5_XX.job
qsub -q i5.q -pe mpich 12 B18_E20_XX.job
qsub -q i6.q -pe mpich 12 B18_E20.5_XX.job
qsub -q i7.q -pe mpich 12 B18_E21_XX.job

qsub -q i8.q -pe mpich 12 B20_E18.job
qsub -q i9.q -pe mpich 12 B20_E18.5.job
