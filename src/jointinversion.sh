#! /bin/sh
#####
# clean up
#####
joint96 39
# define damping
#####
joint96 32 0.5
#####
# Select differential smoothing
#####
joint96 36 1
#####
# set joint weighting between the two data sets
#####
joint96 43 1
#####
# set up repeated run for 5 iterations
#####
joint96 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2
#####
# save last model
#####
joint96 28 modl.out
#####
# plot the model
#####
srfphv96
plotnps -EPS -K -F7 -W10 < SRFPHV96.PLT > figjnt1.eps
