from __future__ import print_function
if __name__ == '__main__':
    print('Loading modules',end='\r')
import modules.glbls as glo
import modules.input as inp
import modules.electronic as el
import modules.vibrational as vi
import modules.output as out
from copy import deepcopy
from time import time

class VHEGEN:
    def __init__(self,problem_dct):
        #initializes VHEGEN input parameters
        self.symmetry = problem_dct['sym']
        self.states = problem_dct['states'] 
        self.modes = problem_dct['modes']
        self.orders = problem_dct['o']
        self.filename = problem_dct['f']
        self.e_coords = 'both' #default e coordinate system
        self.basis = 'both' #default E state basis

    def set_basis(self,basis):
        if 'E' in [i[0] for i in self.states]:
            if basis in ['complex','real','both']:
                self.basis = basis 
            else:
                raise Exception('EBasisError: E basis input must be "real", "complex", or "both".')

    def set_e_coordinates(self,coord_system): #can be 'pol' or 'cart' 
        if hasattr(self,"expansions"):
            raise Exception('ECoordError: e mode coordinate system must be specified before expansion.')
        else:
            if 'E' in [i[0] for i in self.modes]:
                if coord_system in ['pol','cart','both']:
                    self.e_coords = coord_system
                else:
                    raise Exception('ECoordError: e mode coordinate system must be "pol", "cart", or "both".')

    def return_init(self):
        #returns formatted input; for logging
        return inp.return_problem(self.symmetry,self.states,self.modes,self.orders,self.filename)

    def get_eigenvals(self):
        #acquires independent matrix elements and their eigenvalues
        self.eigenvals = el.get_eigenvals(self.symmetry, self.states)

    def return_eigenvals(self):
        #returns formatted eigenvals; for logging
        string='Independent matrix elements: Eigenvals\n\n'
        for e in self.eigenvals:
            string += ('H_'+str(e)+': '+str(self.eigenvals[e]))+'\n'
        return string

    def get_constraints(self):
        self.constraints = vi.load_constraints(self.symmetry,self.eigenvals,self.modes)

    def get_formulas(self):
        self.__formulas = {}
        self.formulas = {}
        for e in self.eigenvals:
            formula = vi.get_root_formula(self.eigenvals[e],self.modes,self.symmetry[1])
            if len(self.modes) < 2:
                formula = vi.adapt_to_unimodal(formula,self.modes)
            self.__formulas[e] = formula
            self.formulas[e] = vi.get_symbolic_formula(formula)
        self.get_constraints()
        self.m_e_info()

    def return_formulas(self):
        string = ('Root formulas:\n\n')
        for e in self.formulas:
            string += ('H_'+str(e)+' : ' + str(self.formulas[e])+'\n')
        return string

    def return_constraints(self):
        return out.format_constraints(self.constraints)

    def get_matrix_form(self):
        self.matrix = el.VibronicMatrix(self.states)
        self.matrix.get_dependencies()

    def get_expansions(self):
        expansions = {}
        count = {}
        parameters = {}
        for o in self.orders:
            print('Order: '+str(o))
            expansions[o] = {}
            count[o] = {}
            parameters[o] = {}
            for e in self.eigenvals: #generate independent matrix elements
                print('Expanding H_'+str(e))
                #for e_coords='both', get_matrix_element_expansions returns normal structure but list of expansions for each matrix element, 1st being pol, 2nd cart.
                expansions[o][e],parameters[o][e] = vi.get_matrix_element_expansion(self.__formulas[e],o,self.constraints[e],self.e_coords,self.eigenvals[e])
                count[o][e] = len(parameters[o][e])
        if hasattr(self,'matrix'):
            expansions = el.get_dependent_elements(expansions,count,self) #map dependent matrix elements

        self.expansions = expansions
        self.parameters = parameters
        self.count = count
        if self.basis == 'real' or self.basis == 'both':
            out.log_append(self.basis_transformation())

    def return_termcount(self):
        if hasattr(self,'count'):
            count_string = ''
            if self.basis == 'complex' or self.basis == 'both':
                count_string += '\nNumber of fitting parameters (complex): \n\n'
                for o in self.count:
                    count_string += 'Order '+str(o)+':\n'
                    for e in self.count[o]:
                        count_string += 'H_'+str(e) +': '+str(self.count[o][e]) +'\n'
            if self.basis == 'real' or self.basis == 'both':
                count_string += '\nNumber of fitting parameters (real): \n\n'
                for o in self.real_count:
                    count_string += 'Order '+str(o)+':\n'
                    for e in self.real_count[o]:
                        count_string += 'H_'+str(e) +': '+str(self.real_count[o][e]) +'\n'
            return count_string
        else:
            raise Exception('TermCountError: Can only count terms after expansion.')
        
    def return_expansions(self):
        two_coords = False
        if len(list(list(self.expansions.values())[0].values())[0]) > 1:
            two_coords = True
        string = 'Matrix element expansions: \n\n'
        for o in self.orders:
            string += 'Order: '+str(o)+'\n\n'
            if two_coords == True:
                string += ('Polar coordinates:\n\n')
            for e in self.expansions[o]:                    
                string += 'H_'+str(e)+' = '+str(self.expansions[o][e][0])+'\n\n'
            if two_coords == True:
                string += ('Cartesian coordinates:\n\n')
                for e in self.expansions[o]:
                    string += 'H_'+str(e)+' = '+str(self.expansions[o][e][1])+'\n\n'
        return string

    def basis_transformation(self):
        if 'E' in [i[0] for i in self.states]:
            if self.basis == 'real' or self.basis =='both':
                if self.basis =='real':
                    mapping = self.matrix.change_basis()
                else:
                    self.matrix_real = deepcopy(self.matrix)
                    mapping = self.matrix_real.change_basis()
                expansions_real,real_count = el.map_elements(self.expansions,mapping)
                self.real_count = real_count
            else:
                return ''
            if self.basis == 'real':
                self.expansions = expansions_real
            elif self.basis == 'both':
                self.expansions_real = expansions_real
            map_string = 'Complex to real basis transformation:\n\n'
            for e in mapping:
                map_string += str(e) +' = ' + str(mapping[e]) +'\n\n'
            return map_string

    def auto(self):
        self.get_eigenvals()
        self.get_matrix_form()
        self.get_formulas()
        self.get_expansions()
        self.basis_transformation()

    def m_e_info(self):
        me_info = {}
        for e in self.eigenvals:
            me_info[e] = [self.eigenvals[e],self.formulas[e],self.constraints[e]]
        self.matrix_element_info = me_info

    def convert(self,sympy_expr,syntax): #valid syntaxs include: 'LaTeX'
        converted_expansions = out.convert_syntax(sympy_expr,syntax)
        return converted_expansions

    def pdflatex(self,path='outputs'):
        TeX_expansions = ''

        for o in self.orders:
            TeX_expansions += R'\subsection{Order: '+str(o)+R'}'+'\n'
            if self.basis == 'real':
                TeX_expansions += R'Number of fitting parameters: ' + out.realcount_format2TeX(self.real_count[o]) + '\n'
            else:
                TeX_expansions += R'Number of fitting parameters: ' + out.count_format2TeX(self.count[o]) + '\n'
            if self.e_coords == 'both':
                TeX_expansions += R'\subsubsection*{Polar e-coordinates:}'+'\n'
                for e in self.expansions[o]:
                    TeX_expansions += '\n'+R'\begin{dmath*}'+'\n'+self.convert(el.format_matrix_element(e),'LaTeX')+'^{('+str(o)+')}='+self.convert(self.expansions[o][e][0],'LaTeX')+R'\\'+'\n'+R'\end{dmath*}' +'\n'
                TeX_expansions += R'\subsubsection*{Cartesian e-coordinates:}' + '\n'
                for e in self.expansions[o]:
                    TeX_expansions += '\n'+R'\begin{dmath*}'+'\n'+self.convert(el.format_matrix_element(e),'LaTeX')+'^{('+str(o)+')}='+self.convert(self.expansions[o][e][1],'LaTeX')+R'\\'+'\n'+R'\end{dmath*}' +'\n'
            else:
                for e in self.expansions[o]:
                    TeX_expansions += '\n'+R'\begin{dmath*}'+'\n'+self.convert(el.format_matrix_element(e),'LaTeX')+'^{('+str(o)+')}='+self.convert(self.expansions[o][e][0],'LaTeX')+R'\\'+'\n'+R'\end{dmath*}' +'\n'
        
        TeX_problem = out.problem_format2TeX(self.states,self.modes)

        if self.basis == 'both': #Create two subsections
            TeX_expansions_real = ''
            for o in self.orders:
                TeX_expansions_real += R'\subsection{Order: '+str(o)+R'}'+'\n'
                TeX_expansions_real += R'Number of fitting parameters: ' + out.realcount_format2TeX(self.real_count[o]) + '\n'
                if self.e_coords == 'both':
                    TeX_expansions_real += R'\subsubsection*{Polar e-coordinates:}'+'\n'
                    for e in self.expansions_real[o]:
                        TeX_expansions_real += '\n'+R'\begin{dmath*}'+'\n'+self.convert(el.format_matrix_element(e),'LaTeX')+'^{('+str(o)+')}='+self.convert(self.expansions_real[o][e][0],'LaTeX')+R'\\'+'\n'+R'\end{dmath*}' +'\n'
                    TeX_expansions_real += R'\subsubsection*{Cartesian e-coordinates:}' + '\n'
                    for e in self.expansions_real[o]:
                        TeX_expansions_real += '\n'+R'\begin{dmath*}'+'\n'+self.convert(el.format_matrix_element(e),'LaTeX')+'^{('+str(o)+')}='+self.convert(self.expansions_real[o][e][1],'LaTeX')+R'\\'+'\n'+R'\end{dmath*}' +'\n'
                else:
                    for e in self.expansions_real[o]:
                        TeX_expansions_real += '\n'+R'\begin{dmath*}'+'\n'+self.convert(el.format_matrix_element(e),'LaTeX')+'^{('+str(o)+')}='+self.convert(self.expansions_real[o][e][0],'LaTeX')+R'\\'+'\n'+R'\end{dmath*}' +'\n'
            TeX_content = out.compose_TeX_both_bases(self.symmetry,self.matrix.format2TeX(),self.matrix_real.format2TeX(),TeX_problem,TeX_expansions,TeX_expansions_real)
        else:
            TeX_content = out.compose_TeX(self.symmetry,self.matrix.format2TeX(),TeX_problem,TeX_expansions)
        out.TeX_write(self.filename,path,TeX_content)
        out.exec_pdflatex(self.filename,path)

def main():
    config, user_input = inp.read_input()
    vhegen = VHEGEN(user_input)
    out.log_append(vhegen.return_init())
    t1 = time()
    vhegen.set_e_coordinates(config[u'e_coords'])
    vhegen.set_basis(config[u'basis'])
    vhegen.get_eigenvals()
    vhegen.get_matrix_form()
    out.log_append(vhegen.return_eigenvals())
    vhegen.get_formulas()
    out.log_append(vhegen.return_formulas())
    out.log_append(vhegen.return_constraints())
    vhegen.get_expansions()
    count_terms = vhegen.return_termcount()
    out.log_append(count_terms)
    out.log_append(vhegen.return_expansions())
    if config[u'pdf_out'] == 'true':
        vhegen.pdflatex()
    if config[u'log_out'] == 'true':
        out.log_write(vhegen.filename,'outputs') 
    print('\nJob complete after '+str(round(time()-t1,3))+' seconds.')
    
if __name__ == '__main__':
    main()
#EOF