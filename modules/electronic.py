import sympy as sym
import sympy.physics.quantum
import numpy as np 
from copy import copy,deepcopy

def kdelta(p,q):
    if p==q:
        return 1
    else:
        return 0

# this doesnt recognize the prime and double prime states --> need to do a string.replace()
def refl_parity(state):
    if ('2' in state):
        return 2
    elif ('1' in state):
        return 1
    else:
        return 0

def inver_parity(state):
    if ('U' in state) or ("''" in state):
        return -1
    elif ('G' in state) or ("'" in state):
        return 1
    else:
        return 0

def get_eigenvals(symmetry, states):
    #this function is analogous to looking up Table 1 in Hickman et al (2018)
    #eigenval order: [Cn_rot,v_refle(Re),v_refle(Im),inver]
    n = int(symmetry[1])
    if len(states)==2:      
        if n % 2 == 1: #odd case
            #interterm cases
            if ('A' in states[0] and 'A' in states[1]):
                #A+A
                eigenvals = {sym.Symbol('A_alphaA_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]} 
            elif ('E' in states[0] and 'A' in states[1]):
                #E+A
                eigenvals = {sym.Symbol('+A'):[0.5*(-1+1j*sym.sqrt(3)),(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]} 
            elif ('E' in states[0] and 'E' in states[1]):
                #E+E
                eigenvals = {sym.Symbol('+_alpha.+_beta'):[1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)],
                             sym.Symbol('+_alpha.-_beta'):[0.5*(-1-1j*sym.sqrt(3)),1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]} 
        
        elif n % 2 == 0: #even case
            #interterm cases
            if ('A' in states[0] and 'A' in states[1]):
                #A+A
                eigenvals = {sym.Symbol('A_alphaA_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]} 
            elif ('B' in states[0] and 'B' in states[1]):
                #B+B
                eigenvals = {sym.Symbol('B_alphaB_beta'):[1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]}
            elif ('A' in states[0] and 'B' in states[1]):
                #A+B
                eigenvals = {sym.Symbol('AB'):[-1,(-1)**(kdelta(refl_parity(states[0]),refl_parity(states[1]))+1),0,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]} 
            elif ('E' in states[0] and 'A' in states[1]):
                #E+A
                eigenvals = {sym.Symbol('+A'):[1j,(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]} 
            elif ('E' in states[0] and 'B' in states[1]):
                #E+B
                eigenvals = {sym.Symbol('+B'):[-1j,(-1)**kdelta(refl_parity(states[1]),2),(-1)**kdelta(refl_parity(states[1]),1),(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]} 
            elif ('E' in states[0] and 'E' in states[1]):
                #E+E
                eigenvals = {sym.Symbol('+_alpha.+_beta'):[1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)],
                             sym.Symbol('+_alpha.-_beta'):[-1,1,-1,(-1)**(kdelta(inver_parity(states[0]),inver_parity(states[1]))+1)]} 
    else:
        if n % 2 == 1: #odd case
            #intraterm cases 
            if 'A' in states[0]:
                #A
                eigenvals = {sym.Symbol('AA'):[1,1,0,1]}
            elif 'E' in states[0]:
                #E
                eigenvals = {sym.Symbol('++'): [1,1,0,1],
                             sym.Symbol('+-'): [0.5*(-1-1j*sym.sqrt(3)),1,-1,1]}
        elif n % 2 == 0: #even case 
            #intraterm cases
            if 'A' in states[0]:
                #A
                eigenvals = {sym.Symbol('AA'):[1,1,0,1]}
            elif 'B' in states[0]:
                #B
                eigenvals = {sym.Symbol('BB'):[1,1,0,1]}
            elif 'E' in states[0]:
                #E
                eigenvals = {sym.Symbol('++'): [1,1,0,1],
                             sym.Symbol('+-'): [-1,1,-1,1]}
    return eigenvals

class VibronicMatrix:
    def __init__(self,states):
        gammastates = []
        Estates = []
        state_components = []
        for s in states:
            if 'E' in s:
                Estates.append(s)
                state_components + ['+','-']
            else:
                gammastates.append(s)
                state_components + []

        if len(gammastates+Estates) == 1: #intraterm cases
            if len(Estates) == 1: #E 
                state_components = ['+','-']
                self.orig = sym.Matrix([[sym.Symbol('++'),sym.Symbol('+-')],
                                        [sym.Symbol('-+'),sym.Symbol('--')]])
                self.symd = sym.Matrix([[sym.Symbol('++'),sym.Symbol('+-')],
                                        [sym.conjugate(sym.Symbol('+-')),sym.Symbol('++')]])
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
                
        self.bras = sym.Matrix([sym.physics.quantum.state.Bra(c) for c in state_components])
        self.kets = sym.transpose(sym.Matrix([i.dual for i in self.bras]))
        self.form = sym.UnevaluatedExpr(self.kets)*sym.UnevaluatedExpr(self.orig)*sym.UnevaluatedExpr(self.bras)

    def __str__(self):
        return str(self.orig)

    def get_dependencies(self):
        independent_set = set()
        for i in self.symd:
            if not isinstance(i,sym.conjugate) and i != 0:
                independent_set.add(i)
        self.dependencies = {}
        for indep_e in independent_set:
            self.dependencies[indep_e] = []
            for c,e in enumerate(self.symd):
                if e == indep_e:
                    self.dependencies[indep_e].append(self.orig[c])
                if sym.conjugate(e) == indep_e:
                    self.dependencies[indep_e].append(sym.conjugate(self.orig[c]))
        return self.dependencies

    def sub_expansions(self,expansions):
        expanded_matrices = {}
        for o in expansions:
            expanded_matrices[o] = copy(self.orig)
        for key in self.dependencies:
            for m_e in self.dependencies[key]:
                #non conj
                if m_e in self.orig:
                    expanded_matrices[o][self.orig.index(m_e)] = expansions[o][key]
                #conj
                elif sym.conjugate(m_e) in self.orig:
                    expanded_matrices[o][self.orig.index(m_e)] = sym.conjugate(expansions[o][key])
        self.expanded_matrices = expanded_matrices 
        return self.expanded_matrices

    def format2TeX(self):
        texmatrix = deepcopy(self.orig)
        for c,i in enumerate(texmatrix):
            if i != 0:
                texmatrix[c] = format_matrix_element(i)
            else:
                texmatrix[c] = i
        form = copy(self.form)
        form = form.subs(self.orig,texmatrix)
        return R'$\hat{H}='+sym.printing.latex(form,mat_delim='(',mat_str='matrix')+'$'

    def change_basis(self): #complex to real basis
        for c,element in enumerate(self.symd):
            if type(element) == sym.conjugate:
                self.symd[c] = sym.re(sym.conjugate(element)) - 1j*sym.im(sym.conjugate(element))
            else:
                self.symd[c] = sym.re(element) + 1j*sym.im(element)
        evolved_matrix = (self.unitary*self.symd*sym.physics.quantum.Dagger(self.unitary))
        for c, element in enumerate(evolved_matrix):
            evolved_matrix[c] = sym.simplify(element)
        self.evolved_matrix = evolved_matrix
        #Change of labels:
        new_matrix_form = copy(self.orig)
        for c,e in enumerate(new_matrix_form):
            if e != 0:
                e_str = str(e)
                e_str = e_str.replace('+','X').replace('-','Y')
                new_matrix_form[c] = sym.Symbol(e_str)
        for c,e in enumerate(self.bras):
            e_str = str(e).replace('<','').replace('|','').replace('+','X').replace('-','Y')
            self.bras = self.bras.subs(e,sym.physics.quantum.state.Bra(e_str))
        for c,e in enumerate(self.kets):
            e_str = str(e).replace('>','').replace('|','').replace('+','X').replace('-','Y')
            self.kets = self.kets.subs(e,sym.physics.quantum.state.Ket(e_str))
        self.orig = new_matrix_form
        self.form = sym.UnevaluatedExpr(self.kets.subs({sym.Symbol('+'):sym.Symbol('X'),sym.Symbol('-'):sym.Symbol('Y')}))*sym.UnevaluatedExpr(new_matrix_form)*sym.UnevaluatedExpr(self.bras.subs({sym.Symbol('+'):sym.Symbol('X'),sym.Symbol('-'):sym.Symbol('Y')}))
        mapping = {}
        for c,e in enumerate(self.orig):
            if e != 0:
                e_str = str(e)
                e_str = e_str.replace('+','X').replace('-','Y')
                mapping[sym.Symbol(e_str)] = self.evolved_matrix[c]
        return mapping
        
def get_dependent_elements(expansion_dct,count,obj):
    print('\nApplying expansion dependencies...\n')
    all_expansions_dct = {}
    dependent_count = {}
    for o in expansion_dct:
        all_expansions_dct[o] = {}
        dependent_count[o] = {}
        for e in expansion_dct[o]:
            for dep_e in obj.matrix.dependencies[e]:
                if isinstance(dep_e,sym.conjugate):
                    all_expansions_dct[o][sym.conjugate(dep_e)] = [sym.conjugate(expansion_dct[o][e][c]) for c,i in enumerate(expansion_dct[o][e])]
                    dependent_count[o][sym.conjugate(dep_e)] = count[o][e]
                else:
                    all_expansions_dct[o][dep_e] = expansion_dct[o][e]
                    dependent_count[o][dep_e] = count[o][e]
    return all_expansions_dct#,dependent_count

def format_matrix_element(matrix_element):
    m_e = str(matrix_element)
    m_e = m_e.replace('.','')
    m_e = m_e.replace('alpha',R'{\alpha}')
    m_e = m_e.replace('beta',R'{\beta}')
    return sym.Symbol('H_{'+m_e+'}')

def prune_dependent_elements(real_matrix_elements):
    indep_m_e = []
    m_e = [str(i) for i in real_matrix_elements]
    #E+A: XA, YA, AX, AY --> XA, YA
    if 'XA' in m_e:
        for i in m_e:
            if i[0] == 'X' or i[0] == 'Y':
                indep_m_e.append(sym.Symbol(i))
    #E+E: all --> Xalpha Xbeta, Xalpha, Ybeta
    elif '.' in m_e[0]:
        m_e_lists = []
        for i in m_e:
            m_e_lists.append(i.split('.'))
        for i in m_e_lists:
            if 'alpha' in i[0] and 'X' in i[0]:
                indep_m_e.append(sym.Symbol('.'.join(i)))
    #E: XX, XY, YX, YY --> XX, XY
    elif 'XX' in m_e:
        for i in m_e:
            if i[0] == 'X':
                indep_m_e.append(sym.Symbol(i))
    else:
        indep_m_e = real_matrix_elements

    return indep_m_e

def map_elements(expansions, mapping):
    coords = 1
    if len(list(list(expansions.values())[0].values())[0]) > 1:
        coords = 2
    mapped_expansions = {}
    term_inheritance = {}
    independent_elements = prune_dependent_elements([i for i in mapping])
    for o in expansions:
        sub_tuples = [(e,expansions[o][e][0]) for e in expansions[o].keys()]
        if coords == 2:
            sub_tuples_2 = [(e,expansions[o][e][1]) for e in expansions[o].keys()]
        mapped_expansions[o] = {}
        term_inheritance[o] = {}
        for k in mapping:
            get_inheritance = False
            if k in independent_elements:
                get_inheritance = True
                inheritance_count = {}
            expanded_terms = []
            if type(mapping[k]) == sym.Add:
                k_args = list(mapping[k].args)
            else:
                k_args = [mapping[k]]
            k_args = (k_args, [i.free_symbols for i in k_args])
            for c,element in enumerate(k_args[1]):
                for c2,item in enumerate(sub_tuples):
                    if item[0].free_symbols == element:
                        if coords == 1:
                            substituted_args = [k_args[0][c].subs(item[0],item[1])]
                        else:
                            substituted_args = [k_args[0][c].subs(item[0],item[1]),k_args[0][c].subs(sub_tuples_2[c2][0],sub_tuples_2[c2][1])]
                        expanded_terms.append(substituted_args)
                        term_arguments = substituted_args[0]

                        if type(term_arguments) == sym.Symbol:
                            term_count = 1
                        elif type(term_arguments) == sym.Mul: #handle numeric prefactor and single terms
                            sum_exists = False
                            for part in term_arguments.args:
                                if type(part) == sym.Add:
                                    sum_exists = True
                                    term_count = len(part.args)
                            if sum_exists == False:
                                term_count = 1
                        else: #sym.Add
                            try:
                                term_count = len(substituted_args[0].args)
                            except AttributeError:
                                term_count = 0
                        if get_inheritance == True:
                            inheritance_count[k_args[1][c].pop()] = term_count

            sum_terms = 0
            for term in expanded_terms:
                sum_terms += term[0]
            if coords == 2:
                sum_terms2 = 0
                for term in expanded_terms:
                    sum_terms2 += term[1]
                sum_terms = [sum_terms,sum_terms2]
            else:
                sum_terms = [sum_terms]
            mapped_expansions[o][k] = sum_terms
            if get_inheritance == True:
                term_inheritance[o][k] = inheritance_count
    return mapped_expansions, term_inheritance

#EOF