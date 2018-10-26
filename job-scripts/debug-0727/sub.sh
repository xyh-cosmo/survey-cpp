#!/bin/bash

qsub -q i0.q -pe mpich 12 E19_b17_beta_10.25.job
qsub -q i3.q -pe mpich 12 E19_b17_beta_15.job

qsub -q i5.q -pe mpich 12 E20_b16_beta_10.25.job
qsub -q i6.q -pe mpich 12 E20_b16_beta_15.job

qsub -q i7.q -pe mpich 12 E20_b17_beta_10.25.job
qsub -q i8.q -pe mpich 12 E20_b17_beta_15.job
