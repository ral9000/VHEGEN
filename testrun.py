#VHEGEN testrun script

import vhegen as VHE

problem = VHE.inp.prepare_input('D3h',"E''","e'+e'",'0,3','testrun')

vhegen_instance = VHE.VHEGEN(problem)

print(vhegen_instance.return_init())

vhegen_instance.set_e_coordinates('both')

vhegen_instance.set_basis('both')

vhegen_instance.auto()

vhegen_instance.pdflatex()

print('\nTestrun complete without errors.')
