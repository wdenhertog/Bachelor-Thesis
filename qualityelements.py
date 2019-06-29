import FEM
from main_programv2 import linear, circle

spfr, efr = FEM.calc_spfr_efr(0.0625, [circle])
m_fem, quality = FEM.fem_meshfit(0.0625, 0.01, 1.722, 1, [circle], spfr)
m_std = FEM.mesh(0.0625, 0.01)
qualityelements = []
efrqualityelements = []
for i in range(len(m_fem.elements)):
    qualityelements.append((1 + ((m_fem.elements[i].size - m_std.elements[i].size) / m_std.elements[i].size) ** 2) * (
            1 + (m_fem.elements[i].skewness - m_std.elements[i].skewness) ** 2) - 1)
for i in range(len(efr.elements)):
    efrqualityelements.append((1 + ((efr.elements[i].size - m_std.elements[i].size) / m_std.elements[i].size) ** 2) * (
            1 + (efr.elements[i].skewness - m_std.elements[i].skewness) ** 2) - 1)
m_fem.plotzeroz('figures_FEM/fem_circle_quality_625.png', mode=3, qualityelements=qualityelements,
                efrqualityelements=efrqualityelements)
efr.plotzeroz('figures_FEM/efr_circle_quality_625.png', mode=3, qualityelements=efrqualityelements,
              efrqualityelements=qualityelements)
