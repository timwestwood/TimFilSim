#!/bin/bash
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=1:mem=6gb
module load tools/prod
module load matlab/R2021a
CODE_LOC=$PBS_O_WORKDIR/post
(cd $CODE_LOC && matlab -nosplash -nodisplay -nojvm -noFigureWindows -r make_wavelength_and_spacing_sweep_efficiency_figure)
