from sympy import Matrix, Symbol, Add, Mul, UnevaluatedExpr, conjugate, re, im, printing, simplify, sqrt
from sympy.physics.quantum import Ket, Bra

mat = Matrix([[1,1],
              [1j,-1j]])*1/sqrt(2)


print(mat)
print(mat.transpose())

state_components = ['+','-']

ket_row = Matrix([Ket(s) for s in state_components]).transpose()

bra_col = Matrix([Bra(s) for s in state_components])

print(ket_row)

print(bra_col)

form =UnevaluatedExpr(ket_row)*UnevaluatedExpr(mat)*UnevaluatedExpr(bra_col)

print(form)