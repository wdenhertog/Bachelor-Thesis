from main_programv2 import *
import xlsxwriter


def calculate_basisfunctions(e):
    bf1 = np.zeros(shape=3)  # linear basis functions in the form: phi_i(x)=a_i_0+a_i_x*x+a_i_y*y
    bf2 = np.zeros(shape=3)  # for every point in the element there is 1 basis function. The basisfunction array has
    bf3 = np.zeros(shape=3)  # the following shape: (a_i_0,a_i_x,a_i_y).
    bf1[1] = 1 / (2 * e.size) * (e.points[1].coordinates[1] - e.points[2].coordinates[1])
    bf1[2] = 1 / (2 * e.size) * (e.points[2].coordinates[0] - e.points[1].coordinates[0])
    bf1[0] = 1 - bf1[1] * e.points[0].coordinates[0] - bf1[2] * e.points[0].coordinates[1]

    bf2[1] = 1 / (2 * e.size) * (e.points[2].coordinates[1] - e.points[0].coordinates[1])
    bf2[2] = 1 / (2 * e.size) * (e.points[0].coordinates[0] - e.points[2].coordinates[0])
    bf2[0] = 1 - bf2[1] * e.points[1].coordinates[0] - bf2[2] * e.points[1].coordinates[1]

    bf3[1] = 1 / (2 * e.size) * (e.points[0].coordinates[1] - e.points[1].coordinates[1])
    bf3[2] = 1 / (2 * e.size) * (e.points[1].coordinates[0] - e.points[0].coordinates[0])
    bf3[0] = 1 - bf3[1] * e.points[2].coordinates[0] - bf3[2] * e.points[2].coordinates[1]
    return np.array([bf1, bf2, bf3])


def calculate_element_matrices(e, lbda, mu):
    # calculation of the element matrices
    bf = calculate_basisfunctions(e)
    e_sxx = np.zeros(shape=(3, 3))
    e_sxy = np.zeros(shape=(3, 3))
    e_syx = np.zeros(shape=(3, 3))
    e_syy = np.zeros(shape=(3, 3))
    for i in range(3):
        for j in range(3):
            e_sxx[i][j] = e.size * ((lbda + 2 * mu) * bf[j][1] * bf[i][1] + mu * bf[j][2] * bf[i][2])
            e_sxy[i][j] = e.size * (lbda * bf[j][2] * bf[i][1] + mu * bf[j][1] * bf[i][2])
            e_syx[i][j] = e.size * (lbda * bf[j][1] * bf[i][2] + mu * bf[j][2] * bf[i][1])
            e_syy[i][j] = e.size * (mu * bf[j][1] * bf[i][1] + (lbda + 2 * mu) * bf[j][2] * bf[i][2])
    return e_sxx, e_sxy, e_syx, e_syy


def stiffness_matrix(grid, lbda, mu):
    # calculation of the big stifness matrix
    # f=(fx_top, fx_left, fx_bottom, fx_right, fy_top, fy_left, fy_bottom, fy_right)
    sxx = np.zeros(shape=(len(grid.points), len(grid.points)))
    sxy = np.zeros(shape=(len(grid.points), len(grid.points)))
    syx = np.zeros(shape=(len(grid.points), len(grid.points)))
    syy = np.zeros(shape=(len(grid.points), len(grid.points)))
    for e in grid.elements:
        e_sxx, e_sxy, e_syx, e_syy = calculate_element_matrices(e, lbda, mu)
        for i in range(3):
            for j in range(3):
                sxx[e.points[i].index][e.points[j].index] += e_sxx[i][j]
                sxy[e.points[i].index][e.points[j].index] += e_sxy[i][j]
                syx[e.points[i].index][e.points[j].index] += e_syx[i][j]
                syy[e.points[i].index][e.points[j].index] += e_syy[i][j]
    s = np.hstack((np.vstack((sxx, syx)), np.vstack((sxy, syy))))

    # TODO: boundary conditions
    # q = np.zeros(shape=(len(grid.points),1))
    return s


def toexcel(s):
    # just a function to check the layout of the matrix quickly
    workbook = xlsxwriter.Workbook('stiffnessmatrix.xlsx')
    worksheet = workbook.add_worksheet()
    for i in range(len(s)):
        for j in range(len(s)):
            worksheet.write(i, j, s[i][j])
    workbook.close()
    return


def checksymmetric(s):
    # check if the (non-)zero places are symmetric
    # (N.B. this does not necessarily mean that the matrix itself is symmetric!)
    for i in range(len(s)):
        for j in range(len(s)):
            if (s[i][j] == 0 and s[j][i] != 0) or (s[i][j] != 0 and s[j][i] == 0):
                return False
    return True


def meshrelaxation(grid, u, lbda, mu):
    s = stiffness_matrix(grid, lbda, mu)

    # restrict the vertical movement of the bottom and top boundaries and the horizontal movement of the left
    # and right boundaries
    for e in grid.boundarypoints:
        if e.iscorner:
            for i in range(2 * len(grid.points)):
                s[e.index][i] = 0
                s[e.index + len(grid.points)][i] = 0
            s[e.index][e.index] = 1
            s[e.index + len(grid.points)][e.index + len(grid.points)] = 1
        elif e.coordinates[0] == 1 or e.coordinates[0] == -1:
            for i in range(2 * len(grid.points)):
                s[e.index][i] = 0
            s[e.index][e.index] = 1
        else:
            for i in range(2 * len(grid.points)):
                s[e.index + len(grid.points)][i] = 0
            s[e.index + len(grid.points)][e.index + len(grid.points)] = 1
    for i in range(len(grid.zpointongrid)):
        for j in range(2 * len(grid.points)):
            s[grid.zpointongrid[i].index][j] = 0
            s[grid.zpointongrid[i].index + len(grid.points)][j] = 0
        s[grid.zpointongrid[i].index][grid.zpointongrid[i].index] = 1
        s[grid.zpointongrid[i].index + len(grid.points)][grid.zpointongrid[i].index + len(grid.points)] = 1
    for i in range(len(u)):
        if u[i] != 0:
            for j in range(2 * len(grid.points)):
                s[i][j] = 0
            s[i][i] = 1
    return np.linalg.solve(s, u)


def qualitycomparison(femgrid, efrgrid):
    # calculate the quality of the mesh generated with FEM (femgrid) and the mesh generated with Euler Forward Relaxation (efrgrid)
    femstats = femgrid.skewness()
    efrstats = efrgrid.skewness()
    femsizes = femgrid.sizes()
    efrsizes = efrgrid.sizes()
    for i in range(len(femsizes)):
        femstats.append(femsizes[i])
        efrstats.append(efrsizes[i])
    femgof = femgrid.gof()
    efrgof = efrgrid.gof()
    for i in range(len(femgof)):
        femstats.append(femgof[i])
        efrstats.append(efrgof[i])
    return femstats, efrstats


def qualitycomparisonplot(femstats, efrstats, lbda, mu, method):
    # create histograms of the quality comparison and save them on the harddrive
    titles = ['Average skewness', 'Maximum skewness', 'Standard deviation of skewnesses', 'Minimum size',
              'Maximum size', 'Standard deviation of sizes', 'Measure of length', 'Measure of area']
    colors = ['xkcd:salmon', 'xkcd:light blue', 'xkcd:green', 'xkcd:yellow', 'xkcd:lilac', 'xkcd:royal blue']
    for i in range(len(femstats)):
        indices = np.arange(2)
        pyplot.bar(indices, [femstats[i][1], efrstats[i][1]], color=colors)
        pyplot.xticks(indices, ['FEM', 'EFR'])
        pyplot.title(titles[i])
        pyplot.savefig(
            'figures_FEM/' + titles[i].replace(' ', '_') + '_l_' + str(lbda) + '_m_' + str(mu) + '_' + method + '.png')
        pyplot.close()
    return


def qualitymeasure(old_mesh, new_mesh):
    # calculates the quality of the new mesh compared to the old mesh:
    # \sum_{elements} ((oldarea(e) - newarea(e)) / oldarea(e) * (oldskewness(e) - newskewness(e)) / oldskewness(e))^2
    quality = 0
    for i in range(len(old_mesh.elements)):
        oldarea = old_mesh.elements[i].size
        if old_mesh.elements[i].skewness <= 0.0000001:
            oldskewness = new_mesh.elements[i].skewness / 10
        else:
            oldskewness = old_mesh.elements[i].skewness
        newarea = new_mesh.elements[i].size
        newskewness = new_mesh.elements[i].skewness
        quality += ((oldarea - newarea) / oldarea * (oldskewness - newskewness) / oldskewness) ** 2
    return quality


def fem_meshfit(gridsize, threshold, lbda, mu, farr):
    m = mesh(gridsize, threshold)
    standardmesh = mesh(gridsize, threshold)
    m.plot()
    s = stiffness_matrix(m, lbda, mu)
    q = np.zeros((2 * len(m.points), 1))
    methodname = 'FEM'
    newtriang = False
    flip = False
    simplepointrelax = False
    simplegridtozeropoint = False
    fixedpointrelaxation = False
    eulerrelaxation = False
    shortestpath = True
    distplot = True
    redistribute = True
    m.initzeroinfo(farr[0])
    m_fit = simplemeshfit(0.3 * gridsize, 'Shortest path fit with redistribution', gridsize, gridsize / 100, farr,
                          newtriang, flip,
                          distplot,
                          simplepointrelax,
                          simplegridtozeropoint, redistribute)
    m_fit.plot()
    u_x_fit = np.zeros(len(m.points))
    u_y_fit = np.zeros(len(m.points))
    for i in range(len(m.points)):
        u_x_fit[i] = (m_fit.points[i].coordinates[0] - m.points[i].coordinates[0])
        u_y_fit[i] = (m_fit.points[i].coordinates[1] - m.points[i].coordinates[1])
    u_xy = np.hstack((u_x_fit, u_y_fit))
    u = meshrelaxation(m, u_xy, lbda, mu)

    for p in m.points:
        p.coordinates[0] += float(u[p.index])
        p.coordinates[1] += float(u[p.index + len(m.points)])
    m.plot()
    m.setupboundaryplist()
    m.finallevelsetinfoupdate()
    for i in m.points:
        m.points[i.index].iszeropath = m_fit.points[i.index].iszeropath
    for i in range(len(m.edges)):
        m.edges[i].islevelset = m_fit.edges[i].islevelset
        m.edges[i].iszeroedge = m_fit.edges[i].iszeroedge
    for e in m.edges:
        e.length = e.calclength()
    for el in m.elements:
        el.size = el.calcsize()
        el.skewness = el.calcskew()

    m.plotzeroz('figures_FEM/FEMfit.png', mode=2)
    m_fit.plotzeroz('figures_FEM/FEMstartsituation.png', mode=2)
    redistribute = False
    eulerrelaxation = True
    m_efr = simplemeshfit(0.3 * hstep, 'Shortest path fit with Euler forward relaxation', hstep, hstep / 100, farr,
                          newtriang, flip, distplot, simplepointrelax, simplegridtozeropoint, redistribute)
    femstats, efrstats = qualitycomparison(m, m_efr)
    qualitycomparisonplot(femstats, efrstats, lbda, mu, farr[0].__name__)
    quality = qualitymeasure(standardmesh, m)
    print('Quality measure = ' + str(quality))
    return m, quality


def recursiveoptimisation(lower_boundary, upper_boundary, gridsize, threshold, farr, tol=1e-5, h=None, c=None, d=None,
                          fc=None, fd=None):
    invphi = (math.sqrt(5) - 1) / 2  # 1/phi
    invphi2 = (3 - math.sqrt(5)) / 2  # 1/phi^2
    (lower_boundary, upper_boundary) = (min(lower_boundary, upper_boundary), max(lower_boundary, upper_boundary))
    if h is None:
        h = upper_boundary - lower_boundary
    if h <= tol:
        return lower_boundary, upper_boundary
    if c is None:
        c = lower_boundary + invphi2 * h
    if d is None:
        d = lower_boundary + invphi * h
    if fc is None:
        mc, fc = fem_meshfit(gridsize, threshold, c, 1, farr)
    if fd is None:
        md, fd = fem_meshfit(gridsize, threshold, d, 1, farr)
    if fc < fd:
        print('lower boundary: ' + str(lower_boundary) + ' upper boundary: ' + str(d))
        return recursiveoptimisation(lower_boundary, d, gridsize, threshold, farr, tol, h * invphi, c=None, fc=None,
                                     d=c, fd=fc)
    else:
        print('lower boundary: ' + str(c) + ' upper boundary: ' + str(upper_boundary))
        return recursiveoptimisation(c, upper_boundary, gridsize, threshold, farr, tol, h * invphi, c=d, fc=fd, d=None,
                                     fd=None)
