from . import matrix, linear_algebra

# << Error >>

ShapeError = matrix.ShapeError
NotSqaureMatrixError = matrix.NotSquareMatrixError
DecompositionError = linear_algebra.DecompositionError
NoInverseMatrixError = linear_algebra.NoInverseMatrixError

# << Class >>

Shape = matrix.Shape
Matrix = matrix.Matrix
Element_Matrix_Arg = linear_algebra.Element_Matrix_Arg
Element_Matrix = linear_algebra.Element_Matrix

# << Common used Matrix >>

zeros = linear_algebra.ones
ones = linear_algebra.ones
nones = linear_algebra.nones
eye = linear_algebra.eye
diag = linear_algebra.diag
iszeros = linear_algebra.iszeros
isones = linear_algebra.isones
isnones = linear_algebra.isnones
iseye = linear_algebra.iseye
isdiag = linear_algebra.isdiag

# << Elementary transformation >>

row_Ele_Trans = linear_algebra.row_Ele_Trans
column_Ele_Trans = linear_algebra.column_Ele_Trans
R_ET = linear_algebra.R_ET
C_ET = linear_algebra.C_ET

# << Gauss-Jordan Simplification >>

gauss = linear_algebra.gauss
jordan = linear_algebra.jordan
rref = linear_algebra.rref
isgauss = linear_algebra.isgauss
isrref = linear_algebra.isrref
solve = linear_algebra.solve

# << LU & LDU & PLU Decompositions >>

LU = linear_algebra.LU
LDU = linear_algebra.LDU
PLU = linear_algebra.PLU

# << Trace >>

trace = linear_algebra.trace

# << Inverse Matrix >>

isinverse = linear_algebra.isinverse

# << Schur Complement >>

schur_complement = linear_algebra.schur_complement
sherman_morrison_woodbury = linear_algebra.sherman_morrison_woodbury
sherman_morrison = linear_algebra.sherman_morrison
schur_C = linear_algebra.schur_C
S_M_W = linear_algebra.S_M_W
S_M = linear_algebra.S_M

# << Determinant >>

det_def = linear_algebra.det_def
det_expand = linear_algebra.det_expand
