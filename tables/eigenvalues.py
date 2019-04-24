from modules.electronic import kdelta, refl_parity, inver_parity
from sympy import sqrt, Symbol

def return_eigenvals(symmetry, states):

	n = int(symmetry[1])

	added_dummy = False
	if len(states) == 1:
		added_dummy = True
		states.append('dummystate')

	if n % 2 == 1: #Trigonal (odd n-axial) problems

		requirements_dct = {0: ['A'],
							1: ['E'],
							2: ['A','A'],
							3: ['E','A'],
							4: ['E','E']}

		eigenvals_dct = {0: {Symbol('A__A'):[1,1,0,1]},

        				 1: {Symbol('+__+'): [1,1,0,1],
                             Symbol('+__-'): [0.5*(-1-1j*sqrt(3)),1,-1,1]},

						 2: {Symbol('A_alpha__A_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
						   
						 3: {Symbol('+__A'):[0.5*(-1+1j*sqrt(3)),(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
						   
						 4: {Symbol('+_alpha__+_beta'):[1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)],
                             Symbol('+_alpha__-_beta'):[0.5*(-1-1j*sqrt(3)),1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]}
        				   
                        }

	else: #Tetragonal (even n-axial) problems

		requirements_dct = {0: ['A'],
							1: ['B'],
							2: ['E'],
							3: ['A','A'],
							4: ['A','B'],
							5: ['E','A'],
							6: ['E','B'],
							7: ['E','E'],
							8: ['B','B']}

		eigenvals_dct = {0: {Symbol('A__A'):[1,1,0,1]},

						 1: {Symbol('B__B'):[1,1,0,1]},

						 2: {Symbol('+__+'): [1,1,0,1],
                             Symbol('+__-'): [-1,1,-1,1]},

						 3: {Symbol('A_alpha__A_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
							
						 4: {Symbol('A__B'):[-1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},

						 5: {Symbol('+__A'):[1j,(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},
						
						 6: {Symbol('+__B'):[-1j,(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]},

						 7: {Symbol('+_alpha__+_beta'):[1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)],
                             Symbol('+_alpha__-_beta'):[-1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]}, 

						 8: {Symbol('B_alpha__B_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]}
						}
						
	if added_dummy == True:
		del states[-1]
		
	problem_reqs = [s[0] for s in states]

	for k in requirements_dct:
		if requirements_dct[k] == problem_reqs:
			return eigenvals_dct[k]

	print('Error: Could not find eigenvalues for states '+str(states)+'.')
	exit()

#EOF