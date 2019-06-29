import FEM
from main_programv2 import linear, dumbbell, circle, star
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

threshold = 0.01
spfrcirc18, efrcirc18 = FEM.calc_spfr_efr(1 / 8, [circle])
spfrlin18, efrlin18 = FEM.calc_spfr_efr(1 / 8, [linear])
spfrdumb18, efrdumb18 = FEM.calc_spfr_efr(1 / 8, [dumbbell])
spfrstar18, efrstar18 = FEM.calc_spfr_efr(1 / 8, [star])

spfrcirc116, efrcirc116 = FEM.calc_spfr_efr(1 / 16, [circle])
spfrlin116, efrlin116 = FEM.calc_spfr_efr(1 / 16, [linear])
spfrdumb116, efrdumb116 = FEM.calc_spfr_efr(1 / 16, [dumbbell])
spfrstar116, efrstar116 = FEM.calc_spfr_efr(1 / 16, [star])

lbdacirc18 = 1.064
lbdalin18 = 3.647
lbdadumb18 = 2.458
lbdastar18 = 1.596

lbdacirc116 = 1.722
lbdalin116 = 4.634
lbdadumb116 = 2.912
lbdastar116 = 2.129

mcirc18, qcirc18 = FEM.fem_meshfit(1 / 8, threshold, lbdacirc18, 1, [circle], spfrcirc18)
mlin18, qlin18 = FEM.fem_meshfit(1 / 8, threshold, lbdalin18, 1, [linear], spfrlin18)
mdumb18, qdumb18 = FEM.fem_meshfit(1 / 8, threshold, lbdadumb18, 1, [dumbbell], spfrdumb18)
mstar18, qstar18 = FEM.fem_meshfit(1 / 8, threshold, lbdastar18, 1, [star], spfrstar18)

mcirc116, qcirc116 = FEM.fem_meshfit(1 / 16, threshold, lbdacirc116, 1, [circle], spfrcirc116)
mlin116, qlin116 = FEM.fem_meshfit(1 / 16, threshold, lbdalin116, 1, [linear], spfrlin116)
mdumb116, qdumb116 = FEM.fem_meshfit(1 / 16, threshold, lbdadumb116, 1, [dumbbell], spfrdumb116)
mstar116, qstar116 = FEM.fem_meshfit(1 / 16, threshold, lbdastar116, 1, [star], spfrstar116)

femstatscirc18, efrstatscirc18 = FEM.qualitycomparison(mcirc18, efrcirc18, [circle])
femstatslin18, efrstatslin18 = FEM.qualitycomparison(mlin18, efrlin18, [linear])
femstatdumb18, efrstatsdumb18 = FEM.qualitycomparison(mdumb18, efrdumb18, [dumbbell])
femstatsstar18, efrstatsstar18 = FEM.qualitycomparison(mstar18, efrstar18, [star])

femstatscirc116, efrstatscirc116 = FEM.qualitycomparison(mcirc116, efrcirc116, [circle])
femstatslin116, efrstatslin116 = FEM.qualitycomparison(mlin116, efrlin116, [linear])
femstatdumb116, efrstatsdumb116 = FEM.qualitycomparison(mdumb116, efrdumb116, [dumbbell])
femstatsstar116, efrstatsstar116 = FEM.qualitycomparison(mstar116, efrstar116, [star])

standardm18 = FEM.mesh(0.125, threshold)
standardstats18 = standardm18.skewness()
standardsizes18 = standardm18.sizes()
for i in range(len(standardsizes18)):
    standardstats18.append(standardsizes18[i])

standardm116 = FEM.mesh(0.0625, threshold)
standardstats116 = standardm116.skewness()
standardsizes116 = standardm116.sizes()
for i in range(len(standardsizes116)):
    standardstats116.append(standardsizes116[i])

titles = ['Average skewness', 'Maximum skewness', 'Standard deviation of skewnesses', 'Minimum size',
          'Maximum size', 'Standard deviation of sizes']
colors = ['xkcd:salmon', 'xkcd:light blue', 'xkcd:green', 'xkcd:light blue', 'xkcd:green', 'xkcd:light blue',
          'xkcd:green', 'xkcd:light blue', 'xkcd:green']
for i in range(len(standardstats18)):
    indices = np.arange(9)
    plt.bar(indices, [standardstats18[i][1], femstatslin18[i][1], efrstatslin18[i][1], femstatscirc18[i][1],
                      efrstatscirc18[i][1], femstatdumb18[i][1], efrstatsdumb18[i][1], femstatsstar18[i][1],
                      efrstatsstar18[i][1]], color=colors)
    plt.xticks(indices,
               ['Standard mesh', 'FEM linear', 'EFR linear', 'FEM circle', 'EFR circle', 'FEM dumbbell', 'EFR dumbbell',
                'FEM star', 'EFR star'], rotation=30, ha='right')
    plt.title(titles[i])
    # plt.show()
    plt.savefig(
        'figures_FEM/' + titles[i].replace(' ', '_') + '125.png')
    plt.close()
