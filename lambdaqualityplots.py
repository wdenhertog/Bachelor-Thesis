import FEM
from main_programv2 import dumbbell, linear, circle, star, cos_times_sin
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

threshold = 0.01

spfrlin18, efrlin18 = FEM.calc_spfr_efr(0.125, [linear])
spfrlin116, efrlin116 = FEM.calc_spfr_efr(0.0625, [linear])
spfrcirc18, efrcirc18 = FEM.calc_spfr_efr(0.125, [circle])
spfrcirc116, efrcirc116 = FEM.calc_spfr_efr(0.0625, [circle])
statscircle18 = []
statscircle116 = []
statslinear18 = []
statslinear116 = []
lbdalist = []
for i in range(40):
    print(i)
    lbda = 0.5 * (i + 1)
    lbdalist.append(lbda)
    m_circ18, quality_circ18 = FEM.fem_meshfit(1 / 8, threshold, lbda, 1, [circle], spfrcirc18)
    m_lin18, quality_lin18 = FEM.fem_meshfit(1 / 8, threshold, lbda, 1, [linear], spfrlin18)
    m_circ116, quality_circ116 = FEM.fem_meshfit(1 / 16, threshold, lbda, 1, [circle], spfrcirc116)
    m_lin116, quality_lin116 = FEM.fem_meshfit(1 / 16, threshold, lbda, 1, [linear], spfrlin116)

    standardstats_circ18 = m_circ18.skewness()
    standardstats_lin18 = m_lin18.skewness()
    standardstats_circ116 = m_circ116.skewness()
    standardstats_lin116 = m_lin116.skewness()

    sizesstats_circ18 = m_circ18.sizes()
    sizesstats_lin18 = m_lin18.sizes()
    sizesstats_circ116 = m_circ116.sizes()
    sizesstats_lin116 = m_lin116.sizes()
    for j in range(len(sizesstats_circ18)):
        standardstats_circ18.append(sizesstats_circ18[j])
        standardstats_lin18.append(sizesstats_lin18[j])
        standardstats_circ116.append(sizesstats_circ116[j])
        standardstats_lin116.append(sizesstats_lin116[j])

    statscircle18.append(standardstats_circ18)
    statslinear18.append(standardstats_lin18)
    statscircle116.append(standardstats_circ116)
    statslinear116.append(standardstats_lin116)
avgskewcirc18 = []
avgskewcirc116 = []
avgskewlin18 = []
avgskewlin116 = []

maxskewcirc18 = []
maxskewcirc116 = []
maxskewlin18 = []
maxskewlin116 = []

stdskewcirc18 = []
stdskewcirc116 = []
stdskewlin18 = []
stdskewlin116 = []

minsizecirc18 = []
minsizecirc116 = []
minsizelin18 = []
minsizelin116 = []

maxsizecirc18 = []
maxsizecirc116 = []
maxsizelin18 = []
maxsizelin116 = []

stdsizecirc18 = []
stdsizecirc116 = []
stdsizelin18 = []
stdsizelin116 = []

for i in range(len(statscircle18)):
    avgskewcirc18.append(statscircle18[i][0][1])
    avgskewlin18.append(statslinear18[i][0][1])
    avgskewcirc116.append(statscircle116[i][0][1])
    avgskewlin116.append(statslinear116[i][0][1])

    maxskewcirc18.append(statscircle18[i][1][1])
    maxskewlin18.append(statslinear18[i][1][1])
    maxskewcirc116.append(statscircle116[i][1][1])
    maxskewlin116.append(statslinear116[i][1][1])

    stdskewcirc18.append(statscircle18[i][2][1])
    stdskewlin18.append(statslinear18[i][2][1])
    stdskewcirc116.append(statscircle116[i][2][1])
    stdskewlin116.append(statslinear116[i][2][1])

    minsizecirc18.append(statscircle18[i][3][1])
    minsizelin18.append(statslinear18[i][3][1])
    minsizecirc116.append(statscircle116[i][3][1])
    minsizelin116.append(statslinear116[i][3][1])

    maxsizecirc18.append(statscircle18[i][4][1])
    maxsizelin18.append(statslinear18[i][4][1])
    maxsizecirc116.append(statscircle116[i][4][1])
    maxsizelin116.append(statslinear116[i][4][1])

    stdsizecirc18.append(statscircle18[i][5][1])
    stdsizelin18.append(statslinear18[i][5][1])
    stdsizecirc116.append(statscircle116[i][5][1])
    stdsizelin116.append(statslinear116[i][5][1])

plt.plot(lbdalist, avgskewcirc18, label='circle, 0.125')
plt.plot(lbdalist, avgskewcirc116, label='circle, 0.0625')
plt.plot(lbdalist, avgskewlin18, label='linear, 0.125')
plt.plot(lbdalist, avgskewlin116, label='linear, 0.0625')
plt.xlabel("\u03BB")
plt.ylabel("Average skewness")
plt.legend(loc='upper right')
plt.savefig('figures_FEM/avgskew_.png')
plt.close()

plt.plot(lbdalist, maxskewcirc18, label='circle, 0.125')
plt.plot(lbdalist, maxskewcirc116, label='circle, 0.0625')
plt.plot(lbdalist, maxskewlin18, label='linear, 0.125')
plt.plot(lbdalist, maxskewlin116, label='linear, 0.0625')
plt.xlabel("\u03BB")
plt.ylabel("Maximum skewness")
plt.legend(loc='upper right')
plt.savefig('figures_FEM/maxskew_.png')
plt.close()

plt.plot(lbdalist, stdskewcirc18, label='circle, 0.125')
plt.plot(lbdalist, stdskewcirc116, label='circle, 0.0625')
plt.plot(lbdalist, stdskewlin18, label='linear, 0.125')
plt.plot(lbdalist, stdskewlin116, label='linear, 0.0625')
plt.xlabel("\u03BB")
plt.ylabel("Standard deviation of skewness")
plt.legend(loc='upper right')
plt.savefig('figures_FEM/stdskew_.png')
plt.close()

plt.plot(lbdalist, minsizecirc18, label='circle, 0.125')
plt.plot(lbdalist, minsizecirc116, label='circle, 0.0625')
plt.plot(lbdalist, minsizelin18, label='linear, 0.125')
plt.plot(lbdalist, minsizelin116, label='linear, 0.0625')
plt.xlabel("\u03BB")
plt.ylabel("Minimum size")
plt.legend(loc='upper right')
plt.savefig('figures_FEM/minsize_.png')
plt.close()

plt.plot(lbdalist, maxsizecirc18, label='circle, 0.125')
plt.plot(lbdalist, maxsizecirc116, label='circle, 0.0625')
plt.plot(lbdalist, maxsizelin18, label='linear, 0.125')
plt.plot(lbdalist, maxsizelin116, label='linear, 0.0625')
plt.xlabel("\u03BB")
plt.ylabel("Maximum size")
plt.legend(loc='upper right')
plt.savefig('figures_FEM/maxsize_.png')
plt.close()

plt.plot(lbdalist, stdsizecirc18, label='circle, 0.125')
plt.plot(lbdalist, stdsizecirc116, label='circle, 0.0625')
plt.plot(lbdalist, stdsizelin18, label='linear, 0.125')
plt.plot(lbdalist, stdsizelin116, label='linear, 0.0625')
plt.xlabel("\u03BB")
plt.ylabel("Standard deviation of size")
plt.legend(loc='upper right')
plt.savefig('figures_FEM/stdsize_.png')
plt.close()
