import FEM
from main_programv2 import dumbbell, linear, circle, star, cos_times_sin

gridsize = 1 / 8
threshold = 0.01
farr = [linear]
tol = 0.25
spfr, efr = FEM.calc_spfr_efr(gridsize, farr)
lower_boundary = 0
upper_boundary = 25
lower_boundary, upper_boundary, lbdalist = FEM.recursiveoptimisation(lower_boundary, upper_boundary, gridsize,
                                                                     threshold, farr, spfr, efr, [], tol)
m, quality = FEM.fem_meshfit(gridsize, threshold, lbdalist[-1][0], 1, farr, spfr)
femstats, efrstats = FEM.qualitycomparison(m, efr, farr)
FEM.qualitycomparisonplot(femstats, efrstats, lbdalist[-1][0], 1, farr[0].__name__, gridsize, threshold, farr)
print('Optimal lambda = ' + str(lbdalist[-1][0]))
