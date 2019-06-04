try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
import matplotlib.pyplot as pyplot
import matplotlib

allmethodnames = ['Cut method', 'Cut and flip method', 'Shortest path fit', 'Shortest path fit with redistribution',
                  'Shortest path fit with Euler forward relaxation']
fnames=False
qualityobjects = []
hstep = 0.125
methodname='Shortest path fit with Euler forward relaxation'
loadname = 'Saved_qualitydata/' + str(round(hstep, 8)) + '/' + methodname.replace(" ", "_") + '.pkl'
try:
    with open(loadname, 'rb') as input:
        qdata = pickle.load(input)
        qualityobjects.append(qdata)
except:
    print('Quality dataset of method <' + methodname + '> is not detected.')
if fnames == False and len(qualityobjects) > 0:
    fnames = qualityobjects[0].fnames
else:
    print('Error, no quality object files found')

nfuns = len(fnames)
nbars = len(qualityobjects[0].after[0]) + 1
width = 0.7 / nbars
colors = ['xkcd:salmon', 'xkcd:light blue', 'xkcd:green', 'xkcd:yellow', 'xkcd:lilac', 'xkcd:royal blue']
plegend = ['Before']

for m, beforedataset in enumerate(qualityobjects[0].after[0]):  # we plot a figure for each measuretype with index m
    measurename = beforedataset[0]  # every datapackage in qualityobjects is a list of the measurename and the data
    beforedata = beforedataset[1]
    fig = pyplot.figure(figsize=(8, 8))
    for f, fun in enumerate(fnames):
        if nbars % 2 == 0:
            plist = [pyplot.bar(f - nbars / 2 * width, beforedata, width, color=colors[0], align='edge')]
            for q, qdata in enumerate(qualityobjects):
                plegend.append(qdata.methodname)
                if fun in qdata.fnames:
                    findex = qdata.fnames.index(fun)
                    plist.append(pyplot.bar(f + (q + 1 - nbars / 2) * width, qdata.after[findex][m][1], width,
                                            color=colors[q + 1], align='edge'))
        else:
            plist = [pyplot.bar(f - nbars // 2 * width, str(beforedata), width, color=colors[0], align='center')]
            for q, qdata in enumerate(qualityobjects):
                plegend.append(qdata.methodname)
                if fun in qdata.fnames:
                    findex = qdata.fnames.index(fun)
                    plist.append(pyplot.bar(f + (q + 1 - nbars // 2) * width, qdata.after[findex][m][1], width,
                                            color=colors[q + 1], align='center'))

    pyplot.title(measurename)
    pyplot.xticks(range(len(qualityobjects[0].fnames)), qualityobjects[0].fnames)
    pyplot.ylabel(measurename)
    xfmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((-3, 4))
    pyplot.gca().yaxis.set_major_formatter(xfmt)
    pyplot.legend(plist, plegend)

    savename = 'Saved_figures/' + str(round(hstep, 8)) + '/Quality_' + measurename.replace(" ", "_") + '.png'
    fig.savefig(savename, dpi=fig.dpi, bbox_inches="tight")
