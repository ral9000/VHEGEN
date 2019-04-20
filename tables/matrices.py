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

			  2: [Matrix([[0,0,'+A'],
			  			  [0,0,'-A'],
			  			  ['A+','A-',0]]),

			  	  Matrix([[0,0,'+A'],
			  			  [0,0,conjugate('+A')],
			  			  [conjugate('+A'),'+A',0]]),

			  	  Matrix([[1,1,0],
			  	  		  [1j,-1j,0],
			  	  		  [0,0,sqrt(2)]])*1/sym.sqrt(2)],

			  3: [Matrix([[0,0,'+B'],
			  			  [0,0,'-B'],
			  			  ['B+','B-',0]]),

			  	  Matrix([[0,0,'+B'],
			  			  [0,0,conjugate('+B')],
			  			  [conjugate('+B'),'+B',0]]),

			  	  Matrix([[1,1,0],
			  	  		  [1j,-1j,0],
			  	  		  [0,0,sqrt(2)]])*1/sym.sqrt(2)],

			  4: [Matrix([[0,0,'+_alpha+_beta','+_alpha-_beta'],
                          [0,0,'-_alpha+_beta','-_alpha-_beta'],
                          ['+_beta+_alpha','+_beta-_alpha',0,0],
                          ['-_beta+_alpha','-_beta-_alpha',0,0]]),

			  	  Matrix([[0,0,'+_alpha+_beta','+_alpha-_beta'],
                          [0,0,conjugate('+_alpha-_beta'),conjugate('+_alpha+_beta')],
                          [conjugate('+_alpha+_beta'),'+_alpha-_beta',0,0],
                          [conjugate('+_alpha-_beta'),'+_alpha+_beta',0,0]]),

			  	  Matrix([[1,1,0,0],
			  	  		  [1j,-1j,0,0],
			  	  		  [0,0,1,1],
			  	  		  [0,0,1j,-1j]])*1/sym.sqrt(2)],

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

			  8: [Matrix([['+,+', '+,-'],
                          ['-,+', '-,-']]),

			  	  Matrix([['++', '+-'],
                          [conjugate('+-'), '++']]),

			  	  Matrix([[1,1],
			  	  		  [1j,-1j]])*1/sym.sqrt(2)]

			  }

def return_matrix(states):
	requirements = [s[0] for s in states]
	for k in requirements_dct:
		if requirements_dct[k] == requirements:
			return matrix_dct[k]
	print('Error: Could not find matrix form for states '+str(states)+'.')
	exit()

