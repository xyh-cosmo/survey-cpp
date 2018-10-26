#!/bin/bash

qsub -q i0.q -pe mpich 12 B18_E20.job  
qsub -q i1.q -pe mpich 12 B20_E18.job
qsub -q i2.q -pe mpich 12 B20_E20.job  
qsub -q i3.q -pe mpich 12 B20_E22.job  
qsub -q i4.q -pe mpich 12 B22_E20.job  
qsub -q i5.q -pe mpich 12 B22_E22.job
