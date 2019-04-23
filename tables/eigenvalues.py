from modules.electronic import kdelta, refl_parity, inver_parity
from sympy import sqrt, Symbol

def return_eigenvals(symmetry, states):

	n = int(symmetry[1])

	if n % 2 == 1: #Trigonal (odd n-axial) problems

		requirements_dct = {0: ['A','A'],
							1: ['E','A'],
							2: ['E','E'],
							3: ['A'],
							4: ['E']}

		eigenvals_dct = {0: {Symbol('A_alphaA_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
						   
						 1: {Symbol('+A'):[0.5*(-1+1j*sqrt(3)),(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
						   
						 2: {Symbol('+_alpha.+_beta'):[1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)],
                             Symbol('+_alpha.-_beta'):[0.5*(-1-1j*sqrt(3)),1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]}, 
        				   
        				 3: {Symbol('AA'):[1,1,0,1]},

        				 4: {Symbol('++'): [1,1,0,1],
                             Symbol('+-'): [0.5*(-1-1j*sqrt(3)),1,-1,1]}
                        }

	else: #Tetragonal (even n-axial) problems

		requirements_dct = {0: ['A','A'],
							1: ['A','B'],
							2: ['E','A'],
							3: ['E','B'],
							4: ['E','E'],
							5: ['B','B'],
							6: ['A'],
							7: ['B'],
							8: ['E']}

		eigenvals_dct = {0: {Symbol('A_alphaA_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
							
						 1: {Symbol('AB'):[-1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},

						 2: {Symbol('+A'):[1j,(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
							
						 3: {Symbol('+B'):[-1j,(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},

						 4: {Symbol('+_alpha.+_beta'):[1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)],
                             Symbol('+_alpha.-_beta'):[-1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]}, 

						 5: {Symbol('B_alphaB_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
							
						 6: {Symbol('AA'):[1,1,0,1]},

						 7: {'BB':[1,1,0,1]},

						 8: {'++': [1,1,0,1],
                             '+-': [-1,1,-1,1]}
						}

	problem_reqs = [s[0] for s in states]

	for k in requirements_dct:
		if requirements_dct[k] == problem_reqs:
			return eigenvals_dct[k]

	print('Error: Could not find eigenvalues for states '+str(states)+'.')
	exit()

#EOF