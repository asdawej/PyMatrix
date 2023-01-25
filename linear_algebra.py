from matrix import *


def zeros(*args: int) -> Matrix:
    '''
    Create a Matrix with 0\n
    ---
    zeros(2, 3) ->\n
    [[0, 0, 0],\n
    [0, 0, 0]]\n
    zeros(2) ->\n
    [[0, 0],\n
    [0, 0]]
    '''
    temp: list[list] = []
    if len(args) == 1:
        for i in range(args[0]):
            temp.append([])
            for _ in range(args[0]):
                temp[-1].append(0)
        return Matrix(temp)
    else:
        for i in range(args[0]):
            temp.append([])
            for _ in range(args[1]):
                temp[-1].append(0)
        return Matrix(temp)


def iszeros(mat: Matrix) -> bool:
    '''
    To check whether it is a zeros Matrix
    '''
    for i in range(1, mat.shape.m+1):
        for j in range(1, mat.shape.n+1):
            if mat[i, j] != 0:
                return False
    return True


def ones(*args: int) -> Matrix:
    '''
    Create a Matrix with 1\n
    ---
    ones(2, 3) ->\n
    [[1, 1, 1],\n
    [1, 1, 1]]\n
    ones(2) ->\n
    [[1, 1],\n
    [1, 1]]
    '''
    if len(args) == 1:
        temp = zeros(args[0])
    else:
        temp = zeros(args[0], args[1])
    temp[...] = 1
    return temp


def isones(mat: Matrix) -> bool:
    '''
    To check whether it is a ones Matrix
    '''
    for i in range(1, mat.shape.m+1):
        for j in range(1, mat.shape.n+1):
            if mat[i, j] != 1:
                return False
    return True


def nones(*args: int) -> Matrix:
    '''
    Create a Matrix with None\n
    ---
    nones(2, 3) ->\n
    [[None, None, None],\n
    [None, None, None]]\n
    nones(2) ->\n
    [[None, None],\n
    [None, None]]
    '''
    if len(args) == 1:
        temp = zeros(args[0])
    else:
        temp = zeros(args[0], args[1])
    temp[...] = None
    return temp


def isnones(mat: Matrix) -> bool:
    '''
    To check whether it is a nones Matrix
    '''
    for i in range(1, mat.shape.m+1):
        for j in range(1, mat.shape.n+1):
            if mat[i, j] != None:
                return False
    return True


def eye(*args: int) -> Matrix:
    '''
    Create a Matrix with only the diagnal being 1\n
    ---
    eye(2, 3) ->\n
    [[1, 0, 0],\n
    [0, 1, 0]]\n
    eye(2) ->\n
    [[1, 0],\n
    [0, 1]]
    '''
    if len(args) == 1:
        temp = zeros(args[0])
    else:
        temp = zeros(args[0], args[1])
    for i in range(1, min(temp.shape.m, temp.shape.n)+1):
        temp[i, i] = 1
    return temp


def iseye(mat: Matrix) -> bool:
    '''
    To check whether it is a eye Matrix
    '''
    for i in range(1, mat.shape.m+1):
        for j in range(1, mat.shape.n+1):
            if i == j:
                if mat[i, j] != 1:
                    return False
            else:
                if mat[i, j] != 0:
                    return False
    return True


class Element_Matrix_Arg(object):
    '''Embedded class, may not use it'''

    def __init__(self, i: int, j: int, k: int) -> NoReturn:
        self.i = i
        self.j = j
        self.k = k

    def __eq__(self, other: 'Element_Matrix_Arg') -> bool:
        return self.i == other.i and self.j == other.j and self.k == other.k


class Element_Matrix(Matrix):
    '''
    Members are /mat/, /shape.m/, /shape.n/, /arg.i/, /arg.j/, /arg.k/\n
    Methods are /T/, /enblock/, /unblock/
    '''

    def __init__(self, n: int, i: int, j: int = None, *, k: Any = None) -> NoReturn:
        '''
        To create an Elementary Matrix, you should give n, i,
        and must give j or k at least one arg\n
        - If i and j, a swap-Matrix for row_i and row_j
        - If i and k, a mul-Matrix for k times row_i
        - If i, j and k, a muladd-Matrix for k times row_i add to row_j
        ---
        Members are /mat/, /shape.m/, /shape.n/, /arg.i/, /arg.j/, /arg.k/
        '''
        temp_pre = eye(n)
        if j:
            if k:
                temp_pre[j, i] = k
            else:
                temp_pre[i, i] = 0
                temp_pre[j, j] = 0
                temp_pre[i, j] = 1
                temp_pre[j, i] = 1
        elif k:
            temp_pre[i, i] = k
        else:
            raise TypeError('Arguments missing for j or k')
        self.mat = temp_pre.mat
        self.shape = Shape(n, n)
        self.arg = Element_Matrix_Arg(i, j, k)


def row_Ele_Trans(
    mat: Matrix,
    i: int, j: int = None, *, k: Any = None,
    record: bool = False
) -> NoReturn | Element_Matrix:
    '''
    Row Elementary Transformation\n
    For i, j, k, please search Elementary_Matrix() for more information
    '''
    if j:
        if k:
            for t in range(1, mat.shape.n+1):
                mat[j, t] += (mat[i, t]*k)
        else:
            for t in range(1, mat.shape.n+1):
                mat[i, t], mat[j, t] = mat[j, t], mat[i, t]
    elif k:
        for t in range(1, mat.shape.n+1):
            mat[i, t] *= k
    else:
        raise TypeError('Arguments missing for j or k')
    if record:
        return Element_Matrix(mat.shape.m, i, j, k=k)


R_ET = row_Ele_Trans


def column_Ele_Trans(
    mat: Matrix,
    i: int, j: int = None, *, k: Any = None,
    record: bool = False
) -> NoReturn | Element_Matrix:
    '''
    Column Elementary Transformation\n
    For i, j, k, please search Elementary_Matrix() for more information
    '''
    if j:
        if k:
            for t in range(1, mat.shape.m+1):
                mat[t, j] += (mat[t, i]*k)
        else:
            for t in range(1, mat.shape.m+1):
                mat[t, i], mat[t, j] = mat[t, j], mat[t, i]
    elif k:
        for t in range(1, mat.shape.m+1):
            mat[t, i] *= k
    else:
        raise TypeError('Arguments missing for j or k')
    if record:
        return Element_Matrix(mat.shape.n, i, j, k=k)


C_ET = column_Ele_Trans


def gauss(
    mat: Matrix,
    record: bool = False,
    preoptimize: bool = True
) -> Matrix | tuple[Matrix, list[Element_Matrix]]:
    '''
    Gauss Simplification, or the partial rref\n
    Pre-optimization is a permutation that\
    put all the maximums of columns on the diagonal\
    to try to simplify calculation\n
    CHANGE and OUTPUT your Matrix\n
    ---
    gauss(
        [[1, 1, 3],\n
        [2, 2, 4]]
    ) ->\n
    [[1, 1, 3],\n
    [0, 0, -2]]\n
    record =
        E(1, 2; -2)
    ->\n
    {preoptimize = True}\n
    [[2, 2, 4],\n
    [0, 0, 1]]\n
    record =
        E(1, 1)\n
        E(1, 2; -1/2)
    '''
    if record:
        path: list[Element_Matrix] = []
        if preoptimize:
            # Pre-optimization
            for i in range(1, min(mat.shape.m, mat.shape.n)+1):
                max_row_save = i
                for t in range(i+1, mat.shape.m+1):
                    if mat[t, i] > mat[max_row_save, i]:
                        max_row_save = t
                if max_row_save != i:
                    path.append(R_ET(mat, i, max_row_save, record=True))
        row_ptr, column_ptr = 1, 1
        while column_ptr <= mat.shape.n:
            # Expected pivot is zero
            if mat[row_ptr, column_ptr] == 0:
                flag = True
                for t in range(row_ptr+1, mat.shape.m+1):
                    # Find a non-zero element, swap the two rows
                    if mat[t, column_ptr] != 0:
                        flag = False
                        path.append(R_ET(mat, row_ptr, t, record=True))
                        break
                # Can not find a non-zero element, move right
                if flag:
                    column_ptr += 1
                    continue
            # Mul-add elimination
            for t in range(row_ptr+1, mat.shape.m+1):
                if mat[t, column_ptr] != 0:
                    path.append(
                        R_ET(
                            mat, row_ptr, t,
                            k=-mat[t, column_ptr]/mat[row_ptr, column_ptr],
                            record=True
                        )
                    )
            row_ptr += 1
            column_ptr += 1
        return mat, path
    else:
        if preoptimize:
            # Pre-optimization
            for i in range(1, min(mat.shape.m, mat.shape.n)+1):
                max_row_save = i
                for t in range(i+1, mat.shape.m+1):
                    if mat[t, i] > mat[max_row_save, i]:
                        max_row_save = t
                if max_row_save != i:
                    R_ET(mat, i, max_row_save)
        row_ptr, column_ptr = 1, 1
        while column_ptr <= mat.shape.n:
            # Expected pivot is zero
            if mat[row_ptr, column_ptr] == 0:
                flag = True
                for t in range(row_ptr+1, mat.shape.m+1):
                    # Find a non-zero element, swap the two rows
                    if mat[t, column_ptr] != 0:
                        flag = False
                        R_ET(mat, row_ptr, t)
                        break
                # Can not find a non-zero element, move right
                if flag:
                    column_ptr += 1
                    continue
            # Mul-add elimination
            for t in range(row_ptr+1, mat.shape.m+1):
                if mat[t, column_ptr] != 0:
                    R_ET(
                        mat, row_ptr, t,
                        k=-mat[t, column_ptr]/mat[row_ptr, column_ptr]
                    )
            row_ptr += 1
            column_ptr += 1
        return mat


def isgauss(mat: Matrix) -> bool:
    '''
    To check whether the Matrix is Gauss-Simplified
    '''
    row_pivots: list[int] = []
    for j in range(1, mat.shape.n+1):
        save = 0
        for i in range(1, mat.shape.m+1):
            if mat[i, j] != 0:
                save = i
        if not save in row_pivots and save != 0:
            row_pivots.append(save)
    t = 1
    for x in row_pivots:
        if x != t:
            return False
        t += 1
    return True


def jordan(
    mat: Matrix,
    record: bool = False
) -> Matrix | tuple[Matrix, list[Element_Matrix]]:
    '''
    Add to Gauss Simplification, the second part of Gauss-Jordan Simplification, or rref\n
    CHANGE and OUTPUT your Matrix\n
    ---
    jordan(
        [[1, 1, 3],\n
        [0, 0, -2]]
    ) ->\n
    [[1, 1, 0],\n
    [0, 0, 1]]\n
    record =
        E(2, 1; 3/2)\n
        E(2; -1/2)
    '''
    if not isgauss(mat):
        raise ValueError('Not row echelon form')
    if record:
        path: list[Element_Matrix] = []
        column_pivots: list[int] = []
        save = 1
        # Find all the pivots
        for j in range(1, mat.shape.n+1):
            if mat[save, j] != 0:
                column_pivots.append(j)
                save += 1
        save -= 1
        # Mul-add elimination
        for x in reversed(column_pivots):
            for i in range(save-1, 0, -1):
                if mat[i, x] != 0:
                    path.append(
                        R_ET(
                            mat, save, i,
                            k=-mat[i, x]/mat[save, x],
                            record=True
                        )
                    )
            save -= 1
        save += 1
        # Convert pivots into ones
        for x in column_pivots:
            if mat[save, x] != 1:
                path.append(
                    R_ET(
                        mat, save,
                        k=1/mat[save, x],
                        record=True
                    )
                )
            save += 1
        return mat, path
    else:
        column_pivots: list[int] = []
        save = 1
        # Find all the pivots
        for j in range(1, mat.shape.n+1):
            if mat[save, j] != 0:
                column_pivots.append(j)
                save += 1
        save -= 1
        # Mul-add elimination
        for x in reversed(column_pivots):
            for i in range(save-1, 0, -1):
                if mat[i, x] != 0:
                    R_ET(
                        mat, save, i,
                        k=-mat[i, x]/mat[save, x]
                    )
            save -= 1
        save += 1
        # Convert pivots into ones
        for x in column_pivots:
            if mat[save, x] != 1:
                R_ET(
                    mat, save,
                    k=1/mat[save, x]
                )
            save += 1
        return mat


def rref(
    mat: Matrix,
    record: bool = False,
    preoptimize: bool = True
) -> Matrix | tuple[Matrix, list[Element_Matrix]]:
    '''
    To get the reduced row echelon form of a Matrix\
    by Gauss-Jordan Simplification\n
    CHANGE and OUTPUT your Matrix\n
    ---
    rref(
        [[1, 1, 3],\n
        [2, 2, 4]]
    ) ->\n
    [[1, 1, 0],\n
    [0, 0, 1]]\n
    record =
        E(1, 2; -2)\n
        E(2, 1; 3/2)\n
        E(2; -1/2)
    ->\n
    {preoptimize = True}\n
    [[1, 1, 0],\n
    [0, 0, 1]]\n
    record =
        E(1, 2)\n
        E(1, 2; -1/2)\n
        E(2, 1; -4)\n
        E(1; 1/2)
    '''
    if record:
        temp, path1 = gauss(mat, True, preoptimize)
        ret, path2 = jordan(temp, True)
        return ret, path1+path2
    else:
        return jordan(gauss(mat, preoptimize=preoptimize))


def isrref(mat: Matrix) -> bool:
    '''
    To check whether the Matrix is rrefed
    '''
    row_pivots: list[int] = []
    for j in range(1, mat.shape.n+1):
        save = 0
        zero_count = 0
        for i in range(1, mat.shape.m+1):
            if mat[i, j] != 0:
                save = i
            else:
                zero_count += 1
        if not save in row_pivots and save != 0:
            if zero_count+1 != mat.shape.m:
                return False
            row_pivots.append(save)
    t = 1
    for x in row_pivots:
        if x != t:
            return False
        t += 1
    return True


if __name__ == '__main__':
    print('zeros & ones & eye test:')
    print(zeros(3), '\n')
    print(ones(2, 3), '\n')
    print(eye(3), '\n')
    print(iszeros(zeros(4, 5)), iseye(eye(7, 2)), '\n')

    print('Element_Matrix test:')
    print(Element_Matrix(4, 1, 4), '\n')
    print(Element_Matrix(4, 1, k=5), '\n')
    print(Element_Matrix(4, 1, 3, k=5), '\n')

    print('Elementary Transformation test:')
    a = Matrix([[1, 2, 3], [4, 5, 6]])
    C_ET(a, 1, 2)
    print(a, '\n')
    C_ET(a, 2, k=3)
    print(a, '\n')
    print(C_ET(a, 1, 2, k=-1, record=True), '\n')

    print('isgauss & isrref test:')
    print(isgauss(Matrix([[1, 2, 0], [0, 3, 4]])))
    print(isgauss(Matrix([[0, 1, 0], [0, 2, 0], [0, 0, 0]])))
    print(isrref(Matrix([[1, 2, 0], [0, 0, 0]])))
    print(isrref(Matrix([[1, 2, 0], [0, 3, 4]])))

    print('gauss test:')
    b = Matrix([[1, 1, 3], [2, 2, 4]])
    c, path = gauss(b[...], True)
    print(c, '\n')
    for i, x in enumerate(path):
        print(i, x, '\n')
    c, path = gauss(b, True, False)
    print(b, '\n')
    for i, x in enumerate(path):
        print(i, x, '\n')

    print('jordan & rref test:')
    b = Matrix([[1, 1, 3], [2, 2, 4]])
    c, path = rref(b[...], True)
    print(c, '\n')
    for i, x in enumerate(path):
        print(i, x, '\n')
    c, path = rref(b, True, False)
    print(b, '\n')
    for i, x in enumerate(path):
        print(i, x, '\n')
