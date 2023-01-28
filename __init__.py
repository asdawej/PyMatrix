from . import matrix, linear_algebra


# << Error >>


ShapeError = matrix.ShapeError


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
iszeros = linear_algebra.iszeros
isones = linear_algebra.isones
isnones = linear_algebra.isnones
iseye = linear_algebra.iseye


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
