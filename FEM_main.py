import FEM
from main_programv2 import dumbbell, linear, circle, star, cos_times_sin

gridsize = 0.125
threshold = 0.01
farr = [circle]
tol = 0.25
lower_boundary = 0
upper_boundary = 20
lower_boundary, upper_boundary = FEM.recursiveoptimisation(lower_boundary, upper_boundary, gridsize, threshold, farr,
                                                           tol)
m, quality = FEM.fem_meshfit(gridsize, threshold, (lower_boundary + upper_boundary) / 2, 1, farr)
