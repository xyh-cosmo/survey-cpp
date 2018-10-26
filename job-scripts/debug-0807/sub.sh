#!/bin/bash

qsub -q x120.q -pe mpich 48 E19_b17_beta_10.25.job-A
qsub -q x120.q -pe mpich 48 E19_b17_beta_10.25.job-B

