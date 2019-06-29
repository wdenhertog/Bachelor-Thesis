import FEM
from main_programv2 import dumbbell, linear, circle, star
import matplotlib.pyplot as plt

threshold = 0.01
spfrcirc18, efrcirc18 = FEM.calc_spfr_efr(1 / 8, [circle])
spfrlin18, efrlin18 = FEM.calc_spfr_efr(1 / 8, [linear])
spfrdumb18, efrdumb18 = FEM.calc_spfr_efr(1 / 8, [dumbbell])
spfrstar18, efrstar18 = FEM.calc_spfr_efr(1 / 8, [star])

spfrcirc116, efrcirc116 = FEM.calc_spfr_efr(1 / 16, [circle])
spfrlin116, efrlin116 = FEM.calc_spfr_efr(1 / 16, [linear])
spfrdumb116, efrdumb116 = FEM.calc_spfr_efr(1 / 16, [dumbbell])
spfrstar116, efrstar116 = FEM.calc_spfr_efr(1 / 16, [star])

qualitycirc18 = []
qualitylin18 = []
qualitydumb18 = []
qualitystar18 = []

qualitycirc116 = []
qualitylin116 = []
qualitydumb116 = []
qualitystar116 = []

lbdalist = []

for i in range(20):
    print(i)
    lbda = i * 0.5
    lbdalist.append(lbda)
    mcirc18, qcirc18 = FEM.fem_meshfit(1 / 8, threshold, lbda, 1, [circle], spfrcirc18)
    mlin18, qlin18 = FEM.fem_meshfit(1 / 8, threshold, lbda, 1, [linear], spfrlin18)
    mdumb18, qdumb18 = FEM.fem_meshfit(1 / 8, threshold, lbda, 1, [dumbbell], spfrdumb18)
    mstar18, qstar18 = FEM.fem_meshfit(1 / 8, threshold, lbda, 1, [star], spfrstar18)

    mcirc116, qcirc116 = FEM.fem_meshfit(1 / 16, threshold, lbda, 1, [circle], spfrcirc116)
    mlin116, qlin116 = FEM.fem_meshfit(1 / 16, threshold, lbda, 1, [linear], spfrlin116)
    mdumb116, qdumb116 = FEM.fem_meshfit(1 / 16, threshold, lbda, 1, [dumbbell], spfrdumb116)
    mstar116, qstar116 = FEM.fem_meshfit(1 / 16, threshold, lbda, 1, [star], spfrstar116)

    qualitycirc18.append(qcirc18)
    qualitylin18.append(qlin18)
    qualitydumb18.append(qdumb18)
    qualitystar18.append(qstar18)

    qualitycirc116.append(qcirc116)
    qualitylin116.append(qlin116)
    qualitydumb116.append(qdumb116)
    qualitystar116.append(qstar116)

plt.plot(lbdalist, qualitycirc18, label='circle, 0.125')
plt.plot(lbdalist, qualitylin18, label='linear, 0.125')
plt.plot(lbdalist, qualitydumb18, label='dumbbell, 0.125')
plt.plot(lbdalist, qualitystar18, label='star, 0.125')
plt.xlabel("\u03BB")
plt.ylabel("Q(\u03BB)")
plt.legend(loc='upper right')
plt.savefig('figures_FEM/totalqualitylambda0.125.png')
plt.close()

plt.plot(lbdalist, qualitycirc116, label='circle, 0.0625')
plt.plot(lbdalist, qualitylin116, label='linear, 0.0625')
plt.plot(lbdalist, qualitydumb116, label='dumbbell, 0.0625')
plt.plot(lbdalist, qualitystar116, label='star, 0.0625')
plt.xlabel("\u03BB")
plt.ylabel("Q(\u03BB)")
plt.legend(loc='upper right')
plt.savefig('figures_FEM/totalqualitylambda0.0625.png')
plt.close()
