import FEM
from main_programv2 import dumbbell, linear, circle, star, cos_times_sin

gridsize = 0.125
threshold = 0.01
farr = [circle]
tol = 0.25
spfr, efr = FEM.calc_spfr_efr(gridsize, farr)
lower_boundary = 0
upper_boundary = 100
lower_boundary, upper_boundary, lbdalist = FEM.recursiveoptimisation(lower_boundary, upper_boundary, gridsize,
                                                                     threshold, farr, spfr, efr, [], tol)
m, quality = FEM.fem_meshfit(gridsize, threshold, lbdalist[-1][0], 1, farr, spfr, efr)
femstats, efrstats = FEM.qualitycomparison(m, efr)
FEM.qualitycomparisonplot(femstats, efrstats, lbdalist[-1][0], 1, farr[0].__name__, gridsize, threshold)
print(lbdalist[-1][0])
