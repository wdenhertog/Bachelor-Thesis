import FEM
from main_programv2 import dumbbell, linear, circle, star, cos_times_sin

gridsize = 1 / 8
threshold = 0.01
farr = [linear]
spfr, efr = FEM.calc_spfr_efr(gridsize, farr)
lbda = 4.634
m, quality = FEM.fem_meshfit(gridsize, threshold, lbda, 1, farr, spfr)
m_stand = FEM.mesh(gridsize, threshold)
print('FEM quality = ' + str(quality))
print('EFR quality = ' + str(FEM.qualitymeasure(m_stand,efr,lbda)))
