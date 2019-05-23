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
