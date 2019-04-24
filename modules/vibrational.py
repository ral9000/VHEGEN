import sympy as sym
import itertools
import os
from modules.glbls import replace_all
from modules.input import reorder
from copy import copy,deepcopy

class Term: #ExpansionTerm
    def __init__(self,prefactor,coeff,funcs):
        self.prefactor = prefactor      #prefactor should be any sympy-digestible argument
        self.coeff = coeff              #coeff should be an instance of Coeff
        self.indices = coeff.indices
        self.funcs = funcs              #funcs should be a list of instances of Mono and/or Trig 
        self.form = [self.prefactor,self.coeff] + self.funcs
    #1  
    def update_indices(self,p):
        #need to update them in coeff and funcs
        for c,i in enumerate(self.indices): #count, index
            self.coeff.indices[c] = self.coeff.indices[c].replace(i,str(p[c]))
            for c2,f in enumerate(self.funcs):
                if isinstance(f,Mono):
                    f.arg = f.arg.replace(i,str(p[c]))
                elif isinstance(f,Trig):
                    for c3,a in enumerate(f.args):
                        f.args[c3] = a.replace(i,str(p[c]))

    def refresh_indices(self):
        self.indices = self.coeff.indices

    def update_form(self):
        self.form = [self.prefactor,self.coeff] + self.funcs

    def get_order(self):
        order = 0
        for f in self.funcs:
            if isinstance(f,Mono):
                order += sym.sympify(f.arg)
        self.order = order

    def zero(self):
        self = 0

    def grab_constraints(self,constraint_dct): #checks global constraints, acquires relevant("local") constraints, returns them in index-position basis.
        applicable_constraints = {}
        for c in constraint_dct:
            if 'nz' in constraint_dct[c]:
                applicable_constraints[c] = constraint_dct[c]
            if c == 'all':
                applicable_constraints[c] = constraint_dct[c]
            else: #summing-index based constraints
                if '&' in c:
                    c_parse = c.split('&')
                    c_index = []
                    for i in self.indices:
                        for t in c_parse:
                            if t in i:
                                c_index.append(self.indices.index(i))
                    if len(c_index) == 2:
                        applicable_constraints['&'.join([str(i) for i in c_index])] = constraint_dct[c]
                else:
                    for i in self.indices:
                        if c in i:
                            if self.indices.index(i) in applicable_constraints:
                                applicable_constraints[self.indices.index(i)].add(constraint_dct[c])
                            else:
                                applicable_constraints[self.indices.index(i)] = set([constraint_dct[c]])
        self.local_constraints = applicable_constraints
        return self.local_constraints

    #2  
    def expand(self,order,eigenvals):
        val_list, native_par_list = gen_index_lists(self.indices,order,eigenvals)
        partitions = []
        expansions = [] #index of the partitions matches index of the summing indices
        n_found = False
        m_found = False
        for c,i in enumerate(self.indices):
            if 'n' in i:
                n_found = True
                n_i = c
            if 'm' in i:
                m_found = True
                m_i = c
        if (n_found == True) and (m_found == True) and (eigenvals[2] == 0) and (eigenvals[0] == 1): #n,m constraint for (e+e) (1,0)
            nm_constraint = True
        else:
            nm_constraint = False

        vals_iterprod = list(itertools.product(*val_list))
        pars_interprod = list(itertools.product(*native_par_list))

        for c,n in enumerate(list(itertools.product(*val_list))):
            meets_real_req = True
            term = deepcopy(self)
            term.update_indices(n)
            term.get_order()
            if (term.order == order) and (n not in partitions):
                meets_req = apply_constraints(self.local_constraints,pars_interprod[c],term,self)
                partitions.append(n)
                if nm_constraint == True:
                    meets_real_req = False
                    # n>0, m anything
                    if n[n_i] > 0:
                        meets_real_req = True
                    # n =0, m>= 0
                    elif (n[n_i] == 0) and (n[m_i] >= 0):
                        meets_real_req = True                       
                if meets_req == True and meets_real_req == True:
                    expansions.append(term)
        self.expansions = expansions
        
    #3
    def compile_expansions(self):   #compile expansion term as sympy expr
        collect_cartesian_terms = None
        for func in self.funcs:
            if isinstance(func,Trig):
                trigfunc = func
                if trigfunc.cartesian == True:
                    if len(trigfunc.args) > 1:
                        collect_cartesian_terms = 2
                    else:
                        collect_cartesian_terms = 1
        free_parameters = []
        compiled_expansions = []
        for term in self.expansions:
            term.coeff.compute()
            term.coeff.compile()
            
            for i,func in enumerate(term.funcs):
                func.compute()
                func.compile()
            compiled_term = 1
            for part in term.form:
                if hasattr(part,'sympy_form') == False:
                    compiled_term = compiled_term*part
                else:
                    compiled_term = compiled_term*part.sympy_form
            #if two e modes:
            if collect_cartesian_terms == 2:
                compiled_term = sym.collect(compiled_term,['x_alpha','y_alpha','x_beta','y_beta'],exact=True)
            #else if 1 e mode:
            elif collect_cartesian_terms == 1:
                compiled_term = sym.collect(compiled_term,['x','y'],exact=True) 
            compiled_expansions.append(compiled_term)

            if term.coeff.sympy_form not in free_parameters and compiled_term != 0:
                free_parameters.append(term.coeff.sympy_form)

        self.expansions = compiled_expansions
        self.parameters = free_parameters

    def compile_formula(self):
        deepcopyself = deepcopy(self)
        sym_formula = ''
        if deepcopyself.prefactor != 1:
            prefactor = str(deepcopyself.prefactor)
            prefactor = prefactor.replace('1j','i')
            sym_formula += ' '+prefactor+'*'
        sym_formula += str(deepcopyself.coeff.compile())
        for func in deepcopyself.funcs:
            func.symcomp()
            sym_formula += func.symb
        return sym_formula
    
    def convert_to_cartesian(self):
        for part in self.funcs:
            part.convert_to_cartesian()

    def convert_to_polar(self):
        for part in self.funcs:
            part.convert_to_polar()

    def adapt_to_unimodal(self,modes):

        if modes[0][0] != 'A': #gamma+a to gamma
            del self.coeff.indices[0]
            self.refresh_indices()
            i = 0
            max_i = len(self.funcs) - 1

            while i <= max_i:
                if isinstance(self.funcs[i],Mono):
                    if self.funcs[i].coord == 'z':
                        del self.funcs[i]
                        i -= 1
                        max_i -= 1

                i += 1

        else: #a+a to a
            for c,i in enumerate(self.indices):
                if 'b' in i:
                    del self.coeff.indices[c]
            self.refresh_indices()

            i = 0
            max_i = len(self.funcs) - 1

            while i <= max_i:
                if isinstance(self.funcs[i],Mono):
                    if self.funcs[i].label == '_beta':
                        del self.funcs[i]
                        i -= 1
                        max_i -= 1
                    else:
                        self.funcs[i].label = ''
                i += 1

        self.update_form()

class Coeff:
    def __init__(self,letterstem,indices):
        self.letter = letterstem[0]
        self.stem = letterstem[1]
        self.indices = indices

    def compute(self):
        for i,index_arg in enumerate(self.indices):
            self.indices[i] = sym.sympify(index_arg)

    def compile(self):
        #compile into sympy symbol
        indexstr = '_{'
        for c,s in enumerate(self.indices):
            indexstr+=str(s)
            if c != (len(self.indices) -1):
                indexstr+=','
            else:
                indexstr += '}'
        self.sympy_form = sym.Symbol(self.letter+'^{'+self.stem+'}'+indexstr,real=True)
        return self.sympy_form


class Mono:
    def __init__(self,coord,label,arg):
        self.coord = coord
        self.label = label
        self.arg = arg
        self.cartesian = False

    def compute(self):
        self.arg = sym.simplify(self.arg)

    def convert_to_cartesian(self):
        if self.coord == 'rho':
            self.coord = 'x,y'
            self.cartesian = True

    def convert_to_polar(self):
        if self.coord == 'x,y':
            self.coord = 'rho'
            self.cartesian = False

    def symcomp(self):
        self.symb = '*'+self.coord+self.label + '**('+ self.arg +')'

    def compile(self):
        if self.cartesian != True:
            self.sympy_form = sym.Symbol(self.coord+self.label,real=True)**(self.arg)
        else:
            expr = (sym.sqrt(sym.Symbol('x'+self.label,real=True)**2 + sym.Symbol('y'+self.label,real=True)**2))**self.arg #(sqrt(x**2 + y**2))
            self.sympy_form = sym.factor(sym.expand(expr)) #sym.factor

class Trig:
    def __init__(self,args):
        self.args = args
        self.cartesian = False

    def compute(self):
        for i,arg in enumerate(self.args):
            self.args[i] = sym.simplify(arg)

    def convert_to_cartesian(self):
        self.cartesian = True

    def convert_to_polar(self):
        self.cartesian = False

class Sin(Trig):
    def __init__(self,args):
        Trig.__init__(self,args)

    def symcomp(self):
        if len(self.args) > 1:
            self.symb = '*sin(('+self.args[0]+')'+'*phi_alpha +('+self.args[1]+')*phi_beta)'
        else:
            self.symb = '*sin(('+self.args[0]+')*phi)'

    def compile(self):
        if len(self.args) > 1:
            if self.cartesian != True:
                self.sympy_form = sym.sin(self.args[0]*sym.Symbol('phi_alpha',real=True)+self.args[1]*sym.Symbol('phi_beta',real=True))
            else:
                expr = sym.sin(self.args[0]*sym.atan2(sym.Symbol('y_alpha',real=True),sym.Symbol('x_alpha',real=True))+self.args[1]*sym.atan2(sym.Symbol('y_beta',real=True),sym.Symbol('x_beta',real=True)))
                self.sympy_form = sym.factor(sym.expand(expr,trig=True)) #sym.factor
        else:
            if self.cartesian != True:
                self.sympy_form = sym.sin(self.args[0]*sym.Symbol('phi',real=True))
            else:
                expr = sym.sin(self.args[0]*sym.atan2(sym.Symbol('y',real=True),sym.Symbol('x',real=True)))
                self.sympy_form = sym.factor(sym.expand(expr,trig=True)) #sym.factor

class Cos(Trig):
    def __init__(self,args):
        Trig.__init__(self,args)

    def symcomp(self):
        if len(self.args) > 1:
            self.symb = '*cos(('+self.args[0]+')'+'*phi_alpha +('+self.args[1]+')*phi_beta)'
        else:
            self.symb = '*cos(('+self.args[0]+')*phi)'
    def compile(self):
        if len(self.args) > 1:
            if self.cartesian != True:
                self.sympy_form = sym.cos(self.args[0]*sym.Symbol('phi_alpha',real=True)+self.args[1]*sym.Symbol('phi_beta',real=True))
            else:
                expr = sym.cos(self.args[0]*sym.atan2(sym.Symbol('y_alpha',real=True),sym.Symbol('x_alpha',real=True))+self.args[1]*sym.atan2(sym.Symbol('y_beta',real=True),sym.Symbol('x_beta',real=True)))
                self.sympy_form = sym.factor(sym.expand(expr,trig=True)) #sym.factor
        else:
            if self.cartesian != True:
                self.sympy_form = sym.cos(self.args[0]*sym.Symbol('phi',real=True))
            else:
                expr = sym.cos(self.args[0]*sym.atan2(sym.Symbol('y',real=True),sym.Symbol('x',real=True)))
                self.sympy_form = sym.factor(sym.expand(expr,trig=True)) #sym.factor

def read_index_arg(expr):
    multiplier = 1
    addend = 0
    if expr == '0':
        return multiplier,addend
    else:
        if expr[0].isdigit() == True:
            multiplier = int(expr[0])
        if '+' in expr:
            addend = int(expr[expr.index('+')+1])
        if '-' in expr:
            addend = int('-'+expr[expr.index('-')+1])
    return multiplier,addend

def pos_ints(max):
    return list(range(0,max+1))

def any_ints(max):
    return list(range(-max,max+1))

def gen_index_attrs(expr,max_val):
    multiplier,addend = read_index_arg(expr)
    if 'm' in expr or 'n' in expr:
        vals = [(multiplier*i + addend) for i in any_ints(max_val)]
        min_val = -max_val
        vals = list(filter(lambda i: i <= max_val, vals))
        vals = list(filter(lambda i: i >= min_val, vals))
    else:
        vals = [(multiplier*i + addend) for i in pos_ints(max_val)]
        vals = list(filter(lambda i: i <= max_val, vals))

    pars = [((i - addend)/multiplier)%2 for i in vals]
    return vals,pars

def requires_signswap(rotational_eigenval):
    rot_eigenval_arg = copy(rotational_eigenval)
    if sym.im(rot_eigenval_arg).is_negative:
        return True
    else:
        return False    

from tables.formulas import requirements_dct, return_formula

def get_root_formula(eigenvals,modes,n_arg):
    n = str(n_arg)
    rotational_eigenval = eigenvals[0]
    refl_Im = eigenvals[2]
    searchmodes = [mode[0] for mode in modes]
    if len(modes) < 2:
        searchmodes.append('A')
        searchmodes = reorder(searchmodes)
    if requires_signswap(rotational_eigenval) == True:
        reqs = [sym.re(rotational_eigenval) - sym.im(rotational_eigenval)*1j,searchmodes]
    else:
        reqs = [rotational_eigenval,searchmodes]
    for k in requirements_dct:
        if requirements_dct[k] == reqs:
            formula = return_formula(n,k)
            if requires_signswap == True:
                for term in formula:
                    if isinstance(term,ET):
                        term.prefactor = sym.conjugate(term.prefactor)          
            if refl_Im == 0:
                max_c = len(formula) - 1
                c = 0
                while c <= max_c:
                    if isinstance(formula[c],Term):
                        if type(formula[c].prefactor) == complex:
                            del formula[c]
                            c -= 1
                            max_c -= 1
                    c += 1
            return(formula)
    raise Exception('VHEGENError: Could not find root formula for rotational eigenvalue '+str(rotational_eigenval)+' and modes '+str(modes)+'.')

def gen_index_lists(indices,order,eigenvals):
    val_list = []
    native_par_list = []
    for i in indices: #regular case
            vals,parities = gen_index_attrs(i,order)
            val_list.append(vals)
            native_par_list.append(parities) #native index parity

    if eigenvals[2] == 0: #take real case
        n_found = False
        m_found = False
        for c,i in enumerate(indices):
            if 'n' in i:
                n_found = True
                n_index = c

            if 'm' in i:
                m_found = True
                m_index = c

        if eigenvals[0] == 1:
            if m_found == True and n_found == False:
                c = 0
                while c <= (len(val_list[m_index]) - 1):
                    if val_list[m_index][c] < 0:
                        del val_list[m_index][c]
                        del native_par_list[m_index][c]
                        c -= 1
                    c += 1
            if n_found == True and m_found == False:
            #(e+a): n >= 0
                c = 0
                while c <= (len(val_list[n_index]) - 1):
                    if val_list[n_index][c] < 0:
                        del val_list[n_index][c]
                        del native_par_list[n_index][c]
                        c -= 1
                    c += 1
        elif eigenvals[0] == -1:
            if m_found == True and n_found == False:
                c = 0
                while c <= (len(val_list[m_index]) - 1):
                    if val_list[m_index][c] < 0:
                        del val_list[m_index][c]
                        del native_par_list[m_index][c]
                        c -= 1
                    c += 1
    return val_list, native_par_list
        
def adapt_to_unimodal(formula,modes):
    adapted_formula = []
    if formula != [0]:
        for term in formula:
            term.adapt_to_unimodal(modes)
            adapted_formula.append(term)
        
    return adapted_formula

def get_sym_props(pointgroup):
    sym_props = {'inver': False,
                 'refl': False,
                 'hrefl':False}
    if ('D' in pointgroup) or ('V' in pointgroup):
        sym_props['refl'] = True
    if 'H' in pointgroup:
        if '4' in pointgroup:
            sym_props['inver'] = True
        elif '3' in pointgroup:
            sym_props['hrefl'] = True
    return sym_props

def get_constraint(eigenvals, vib_modes, operation):
    princ_rot = eigenvals[0]
    refl_Re, refl_Im = eigenvals[1], eigenvals[2]
    inver_eigenval = eigenvals[3]
    try:
        os.chdir(os.getcwd() + '/constraints')
        constraintfiles = os.listdir(os.getcwd())
        constraintargs = [sym.sympify(i[0:i.index('_')].replace('X','*')) for i in constraintfiles]
        #Open refl.sym file
        for arg in constraintargs:
            if sym.simplify(princ_rot-arg) == 0:
                with open(replace_all(str(arg),{'*':'X',' ':''})+"_"+operation+".sym","r") as sym_file:
                    constraint_lines  = sym_file.readlines()
            elif sym.simplify(sym.re(princ_rot)-sym.im(princ_rot)*1j - arg) == 0:
                with open(replace_all(str(arg),{'*':'X',' ':''})+"_"+operation+".sym","r") as sym_file:
                    constraint_lines  = sym_file.readlines()
        os.chdir('..')
    except OSError as e:
        print(e)
        
    vib_modes = [i.replace("''",'"') for i in vib_modes]

    if operation == 'refl':
        eigen_req = '[' + str(refl_Re) + ',' + str(refl_Im) + ']'
        vib_modes = [replace_all(i,{'G':'','U':'',"'":'','"':''}) for i in vib_modes]

    elif operation == 'inver' or operation == 'hrefl':
        eigen_req = '[' + str(inver_eigenval) + ']'
        vib_modes = [replace_all(i,{'1':'','2':''}) for i in vib_modes]

    if len(vib_modes) == 1:        
        if operation == 'refl':
            vib_modes.append('A1')
        elif operation == 'inver':
            vib_modes.append('AG')
        elif operation == 'hrefl':
            vib_modes.append("A'")
    mode_req = vib_modes
    max_count = len(constraint_lines) - 1
    for count,line in enumerate(constraint_lines):
        modes_in_line = ''.join(line.split(',')[0])
        modes_in_line = modes_in_line[1:-1]
        modes_in_line = modes_in_line.split('.')
        line = line.replace('\n','')  
        if set([str(mode) for mode in mode_req]) == set(modes_in_line) and eigen_req in line:
            match_line = line
        if (count == max_count):
            try:
                match_line
            except NameError: #no matching constraints
                return {}
    constraints = match_line.split(': ',1)[1]

    constraints = constraints.split(',')
    constraints_dict = {}
    for constraint in constraints:
        index = constraint.split(' ')[0]
        restriction = ' '.join(constraint.split(' ')[1:])
        constraints_dict[index] = restriction
    if 'all' in constraints_dict:
        if constraints_dict['all'] == 'nr':
            del constraints_dict['all'] #check if any composite cases where 
    return constraints_dict
    
def load_constraints(sym, eigenvals, vib_modes):
    sym_props = get_sym_props(sym)
    matrix_element_constraints = {}
    for e in eigenvals:
        constraints_dict = {}
        if sym_props['refl'] == True:
            constraints_refl = get_constraint(eigenvals[e], vib_modes, 'refl')
            for key in constraints_refl:
                if key not in constraints_dict:
                    constraints_dict[key] = constraints_refl[key]
                else:
                    constraints_dict[key] += constraints_refl[key]
        if sym_props['inver'] == True:
            constraints_inver = get_constraint(eigenvals[e], vib_modes, 'inver')
            for key in constraints_inver:
                if key not in constraints_dict:
                    constraints_dict[key] = constraints_inver[key]
                else:
                    constraints_dict[key] += constraints_inver[key]
        if sym_props['hrefl'] == True:
            constraints_hrefl = get_constraint(eigenvals[e], vib_modes, 'hrefl')
            for key in constraints_hrefl:
                if key not in constraints_dict:
                    constraints_dict[key] = constraints_hrefl[key]
                else:
                    constraints_dict[key] += constraints_hrefl[key]
        matrix_element_constraints[e] = constraints_dict
    return matrix_element_constraints

def apply_constraints(local_constraints, p, fitted_term,unfitted_term): #returns whether constraints are sat (True) or not sat (False)
    sat_constraints = 0
    if local_constraints != {}:
        if 'all' in local_constraints.keys():
            if 'na' in local_constraints['all']: #all na case
                return False #return not sat
            else:
                if len(local_constraints) == 1: #all nr case. len 1 check ensures theres no other constraints to consider.
                    return True #return sat

        for c in local_constraints:
            if type(c) == int: #handle summing-index based constraints
                if type(local_constraints[c]) == set:
                    len_cond = len(local_constraints[c])
                    met_cond = 0
                    for i in local_constraints[c]:
                        if constraint_funcs[i](c,p) == True:
                            met_cond += 1
                    if met_cond == len_cond:
                        sat_constraints += 1
                elif constraint_funcs[local_constraints[c]](c,p) == True:
                    sat_constraints += 1
            elif type(c.split('&')[0]) == int: #handle pairwise summing-index based constraints
                if constraint_funcs[local_constraints[c]]([int(i) for i in c.split('&')],p) == True:
                    sat_constraints += 1
            elif 'nz' in local_constraints[c]: #handle nz constraints
                if 'if' in local_constraints[c]: #has condition
                    cond = (local_constraints[c].split(' '))[-2:]
                else:
                    cond = []
                if '&' in c: #pairwise nz
                    c1 = c.split('&')[0]
                    c2 = c.split('&')[1]
                    if nz(c1,cond,fitted_term,unfitted_term,p) == True or nz(c2,cond,fitted_term,unfitted_term,p) == True:
                        sat_constraints +=1
                else:
                    if nz(c,cond,fitted_term,unfitted_term,p) == True:
                        sat_constraints += 1
    return sat_constraints == len(local_constraints)

#Constraint functions
def even(i,p):
    if p[i] % 2 != 0:
        return False
    else:
        return True

def odd(i,p):
    if p[i] % 2 == 0:
        return False
    else:
        return True

def ee_or_oo(i,p):
    if (p[i[0]] % 2 + p[i[1]] % 2) == 1:
        return False
    else:
        return True

def eo_or_oe(i,p):
    if (p[i[0]] % 2 + p[i[1]] % 2) != 1:
        return False
    else:
        return True

def non_neg(i,p):
    if p[i] >= 0:
        return True
    else:
        return False

def nm_postproc(i,p):
    if p[i][0] >  0:
        return True 
    if p[i][0] == 0 and p[i][1] >= 0:
        return True
    else:
        return False 

def na(i,p):
    return False

def nr(i,p):
    return True

def nz(arg,cond,fitted_term,unfitted_term,p): #coeff nz and trig nz constraints
    meets_req = False
    satisfy_arg = False
    if cond == []:
        meets_cond = True
    else: #if condition isn't true, constraint need not be applied.
        for c,i in enumerate(unfitted_term.indices):
            if cond[0] in i:
                cond_index = c
        if cond[1] == 'even':
            if p[cond_index] % 2 == 0:
                meets_cond = True
            else:
                meets_cond = False 
        else: #odd
            if p[cond_index] % 2 != 0:
                meets_cond = True
            else:
                meets_cond = False
    if meets_cond == False: #bypass
        meets_req = True
        return meets_req

    else: #nz constrain
        if arg in ['sin','cos']: #handle trig nz
            for part in fitted_term.funcs:
                if (arg == 'sin') and (isinstance(part,Sin)):
                    satisfy_arg = True
                elif (arg == 'cos') and (isinstance(part,Cos)):
                    satisfy_arg = True
        else: #handle coeff nz 
            satisfy_arg = coeff_nz(arg,fitted_term)
    if satisfy_arg == True:
        meets_req = True
    else:
        meets_req = False
    return meets_req

def coeff_nz(arg,term):
    validate = 0
    if arg[1] != '%':
        if term.coeff.stem == arg[1]:
            validate += 1
    else:
        validate += 1
    if not arg.endswith(R'%%'):
        index_parities = [arg[2],arg[3]]
        for c,i in enumerate(index_parities):
            if i == 'e':
                if int(term.indices[c]) % 2 == 0:
                    validate += 1
            elif i == 'o':
                if int(term.indices[c]) % 2 != 0:
                    validate += 1
    else:
        validate += 2
    if validate == 3:
        return True
    else:
        return False

def get_symbolic_formula(formula):
    symb_formula = ''
    if formula != [0]:
        for c,term_formula in enumerate(formula):
            symb_term_formula = term_formula.compile_formula()
            if symb_term_formula[1] != '-' and c != 0:
                symb_formula += '+'
            symb_formula += term_formula.compile_formula()
    return symb_formula

def get_matrix_element_expansion(formula,order,constraints,e_coord_system,eigenvals):
    expansions = []
    params = []
    if e_coord_system == 'both':
        iterates = 2
    else:
        iterates = 1
    for i in range(iterates):
        total_expansion = []
        if formula != [0]:
            for term_formula in formula:
                term_formula.grab_constraints(constraints)
                if e_coord_system == 'pol' or i == 0:
                    term_formula.convert_to_polar()
                if e_coord_system == 'cart' or i == 1:
                    term_formula.convert_to_cartesian()
                term_formula.expand(order,eigenvals)
                term_formula.compile_expansions()
                total_expansion += term_formula.expansions
                params += term_formula.parameters
        expansions.append(total_expansion)
    params = set(params)
    built_expansions = [0]*iterates
    for c,i in enumerate(expansions):
        for exp in i:
            built_expansions[c] += exp
            
    return built_expansions,params

constraint_funcs = {'even':even,
                    'odd':odd,
                    'ee_or_oo':ee_or_oo,
                    'eo_or_oe':eo_or_oe,
                    'non_neg': non_neg,
                    'nm_postproc': nm_postproc,
                    'na':na,
                    'nr':nr}

#EOF