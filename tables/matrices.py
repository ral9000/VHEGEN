from sympy import Matrix, Symbol, conjugate, sqrt

requirements_dct = {0: ['A','A'],
                    1: ['A','B'],
                    2: ['E','A'],
                    3: ['E','B'],
                    4: ['E','E'],
                    5: ['B','B'],
                    6: ['A'],
                    7: ['B'],
                    8: ['E']}

#Structure of matrix_dct values: [original H ,symmetrized H, transformation U]
matrix_dct = {0: [Matrix([[0,'A_alphaA_beta'],
                          ['A_betaA_alpha',0]]),

                  Matrix([[0,'A_alphaA_beta'],
                          ['A_alphaA_beta',0]]),

                  Matrix([[1,0],
                          [0,1]])],

              1: [Matrix([[0,'AB'],
                          ['BA',0]]),

                  Matrix([[0,'AB'],
                          ['AB',0]]),

                  Matrix([[1,0],
                          [0,1]])],

              2: [Matrix([[0,0,Symbol('+A')],
                          [0,0,Symbol('-A')],
                          [Symbol('A+'),Symbol('A-'),0]]),

                  Matrix([[0,0,Symbol('+A')],
                          [0,0,conjugate(Symbol('+A'))],
                          [conjugate(Symbol('+A')),Symbol('+A'),0]]),

                  Matrix([[1,1,0],
                          [1j,-1j,0],
                          [0,0,sqrt(2)]])*1/sqrt(2)],

              3: [Matrix([[0,0,Symbol('+B')],
                          [0,0,Symbol('-B')],
                          [Symbol('B+'),Symbol('B-'),0]]),

                  Matrix([[0,0,Symbol('+B')],
                          [0,0,conjugate(Symbol('+B'))],
                          [conjugate(Symbol('+B')),Symbol('+B'),0]]),

                  Matrix([[1,1,0],
                          [1j,-1j,0],
                          [0,0,sqrt(2)]])*1/sqrt(2)],

              4: [Matrix([[0,0,Symbol('+_alpha.+_beta'),Symbol('+_alpha.-_beta')],
                          [0,0,Symbol('-_alpha.+_beta'),Symbol('-_alpha.-_beta')],
                          [Symbol('+_beta.+_alpha'),Symbol('+_beta.-_alpha'),0,0],
                          [Symbol('-_beta.+_alpha'),Symbol('-_beta.-_alpha'),0,0]]),

                  Matrix([[0,0,Symbol('+_alpha.+_beta'),Symbol('+_alpha.-_beta')],
                          [0,0,conjugate(Symbol('+_alpha.-_beta')),conjugate(Symbol('+_alpha.+_beta'))],
                          [conjugate(Symbol('+_alpha.+_beta')),Symbol('+_alpha.-_beta'),0,0],
                          [conjugate(Symbol('+_alpha.-_beta')),Symbol('+_alpha.+_beta'),0,0]]),

                  Matrix([[1,1,0,0],
                          [1j,-1j,0,0],
                          [0,0,1,1],
                          [0,0,1j,-1j]])*1/sqrt(2)],

              5: [[Matrix([[0,'B_alphaB_beta'],
                          ['B_betaB_alpha',0]]),

                  Matrix([[0,'B_alphaB_beta'],
                          ['B_alphaB_beta',0]]),

                  Matrix([[1,0],
                          [0,1]])]],

              6: [Matrix([['AA']]),

                  Matrix([['AA']]),

                  Matrix([[1]])],

              7: [Matrix([['BB']]),

                  Matrix([['BB']]),

                  Matrix([[1]])],

              8: [Matrix([[Symbol('++'), Symbol('+-')],
                          [Symbol('-+'), Symbol('--')]]),

                  Matrix([[Symbol('++'), Symbol('+-')],
                          [conjugate(Symbol('+-')), Symbol('++')]]),

                  Matrix([[1,1],
                          [1j,-1j]])*1/sqrt(2)]

              }

def return_matrices(states):
    requirements = [s[0] for s in states]
    for k in requirements_dct:
        if requirements_dct[k] == requirements:
            return matrix_dct[k][0], matrix_dct[k][1], matrix_dct[k][2]
    print('Error: Could not find matrix form for states '+str(states)+'.')
    exit()

