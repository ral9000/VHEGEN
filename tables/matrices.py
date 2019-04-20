from sympy import Matrix, Symbol, conjugate

requirements_dct = {0: ['A','A'],
							1: ['A','B'],
							2: ['E','A'],
							3: ['E','B'],
							4: ['E','E'],
							5: ['B','B'],
							6: ['A'],
							7: ['B'],
							8: ['E']}

#Structure of values: [original H ,symmetrized H, transformation U]
matrix_dct = {0: [],

			  1: [],

			  2: [],

			  3: [],

			  4: [],

			  5: [],

			  6: [],

			  7: [],

			  8: []

}



if len(gammastates+Estates) == 1: #intraterm cases
            if len(Estates) == 1: #E 

                state_components = ['+','-']

                self.orig = sym.Matrix([['+,+', '+,-'],
                                        ['-,+', '-,-']])

                self.symd = sym.Matrix([['+,+', '+,-'],
                                        ['+,-*', '+,+']])

                self.unitary = 1/sym.sqrt(2)*sym.Matrix([[1,1],[1j,-1j]])
            elif len(gammastates) == 1: #gamma
                state_components = [gammastates[0][0]]
                self.orig = sym.Matrix([[sym.Symbol(gammastates[0][0]+gammastates[0][0])]])
                self.symd = self.orig
                self.unitary = sym.Matrix([[1]])

        else: #interterm cases
            if len(Estates) == 2: #E+E
                state_components = ['+_alpha','-_alpha','+_beta','-_beta']
                self.orig = sym.Matrix([[0,0,sym.Symbol('+_alpha.+_beta'),sym.Symbol('+_alpha.-_beta')],
                                        [0,0,sym.Symbol('-_alpha.+_beta'),sym.Symbol('-_alpha.-_beta')],
                                        [sym.Symbol('+_beta.+_alpha'),sym.Symbol('+_beta.-_alpha'),0,0],
                                        [sym.Symbol('-_beta.+_alpha'),sym.Symbol('-_beta.-_alpha'),0,0]])
                self.symd = sym.Matrix([[0,0,sym.Symbol('+_alpha.+_beta'),sym.Symbol('+_alpha.-_beta')],
                                        [0,0,sym.conjugate(sym.Symbol('+_alpha.-_beta')),sym.conjugate(sym.Symbol('+_alpha.+_beta'))],
                                        [sym.conjugate(sym.Symbol('+_alpha.+_beta')),sym.Symbol('+_alpha.-_beta'),0,0],
                                        [sym.conjugate(sym.Symbol('+_alpha.-_beta')),sym.Symbol('+_alpha.+_beta'),0,0]])
                self.unitary = 1/sym.sqrt(2)*sym.Matrix([[1,1,0,0],[1j,-1j,0,0],[0,0,1,1],[0,0,1j,-1j]])

            elif len(Estates) == 1 and len(gammastates) == 1: #E+gamma
                state_components = ['+','-',gammastates[0][0]]
                self.orig = sym.Matrix([[0,0,sym.Symbol('+'+gammastates[0][0])],
                                        [0,0,sym.Symbol('-'+gammastates[0][0])],
                                        [sym.Symbol(gammastates[0][0]+'+'),sym.Symbol(gammastates[0][0]+'-'),0]])

                self.symd = sym.Matrix([[0,0,sym.Symbol('+'+gammastates[0][0])],
                                        [0,0,sym.conjugate(sym.Symbol('+'+gammastates[0][0]))],
                                        [sym.conjugate(sym.Symbol('+'+gammastates[0][0])),sym.Symbol('+'+gammastates[0][0]),0]])
                self.unitary = 1/sym.sqrt(2)*sym.Matrix([[1,1,0],[1j,-1j,0],[0,0,sym.sqrt(2)]])

            elif len(gammastates) == 2: #gamma+gamma
                if gammastates[0][0] == gammastates[1][0]: #need alpha/beta labelling
                    state_components = [gammastates[0][0]+'_alpha',gammastates[1][0]+'_beta']
                    self.orig = sym.Matrix([[0,sym.Symbol(gammastates[0][0]+'_alpha'+gammastates[1][0]+'_beta')],
                                            [sym.Symbol(gammastates[1][0]+'_beta'+gammastates[0][0]+'_alpha'),0]])
                    self.symd = sym.Matrix([[0,sym.Symbol(gammastates[0][0]+'_alpha'+gammastates[1][0]+'_beta')],
                                            [sym.Symbol(gammastates[0][0]+'_alpha'+gammastates[1][0]+'_beta'),0]])
                
                else:
                    state_components = [gammastates[0][0],gammastates[1][0]]
                    self.orig = sym.Matrix([[0,sym.Symbol(gammastates[0][0]+gammastates[1][0])],
                                            [sym.Symbol(gammastates[1][0]+gammastates[0][0]),0]])
                    self.symd = sym.Matrix([[0,sym.Symbol(gammastates[0][0]+gammastates[1][0])],
                                            [sym.Symbol(gammastates[0][0]+gammastates[1][0]),0]])
                self.unitary = sym.Matrix([[1,0],[0,1]])
                