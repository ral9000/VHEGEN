import sympy as sym
import os
import subprocess
import modules.electronic as el

log_str = ''

def log_append(string):
    #add contents to log file
    global log_str
    log_str += str(string) + '\n'

def format_mulliken_scripts(irrep,t):
    if t == 's':
        irrep_str = irrep[0].upper()
    else:
        irrep_str = irrep[0].lower()
    primes = 0
    for char in irrep:
        if char == "'":
            primes += 1
    if primes == 0:
        irrep_str += '_{' + irrep[1:].lower() + '}'
    else:
        irrep_str += primes*"'" + '_{' + irrep[1:].replace("'","") + '}'
    return irrep_str

def problem_format2TeX(states,modes):
    problem_str = '$'
    if len(states) > 1:
        problem_str += '('
    for c,s in enumerate(states):
        if c>0:
            problem_str +='+'
        problem_str += format_mulliken_scripts(s,'s')
    if len(states) > 1:
        problem_str += ')'
    problem_str += R' \otimes '
    if len(modes) > 1:
        problem_str += '('
    for c2,m in enumerate(modes):
        if c2>0:
            problem_str+='+'
        problem_str += format_mulliken_scripts(m,'m')
    if len(modes) > 1:
        problem_str += ')'
    problem_str += '$'
    return problem_str

def prune_dependent_terms(count):
  keys = list(count.keys())
  unique_keys = []
  unique_count = {}
  for i in keys:
    if sym.Symbol(str(i)[::-1]) not in unique_keys:
      unique_keys.append(i)
  for k in unique_keys:
    unique_count[k] = count[k]
  return(unique_count)

def count_format2TeX(count): 
  #{+A: 8, -A: 8, A+: 8, A-: 8}
  count = prune_dependent_terms(count)
  count_str = ''
  for c,m_e in enumerate(count):
    count_str += '$'+str(el.format_matrix_element(m_e))+'$: ' + '$'+str(count[m_e]) +'$'
    if c != len(count)-1:
      count_str += ', '
    else:
      count_str += '.'
  return count_str

def realcount_format2TeX(count): #includes inheritance from complex
  #{XA: {+A: 6, ++: 10}, YA: {+A: 6}, AX: {+A: 6}, AY: {+A: 6}}
  count = prune_dependent_terms(count)
  count_str = ''
  for c,m_e in enumerate(count):
    sumval = sum(count[m_e][i] for i in count[m_e])
    count_str += '$'+str(el.format_matrix_element(m_e))+'$: ' + '$'+str(sumval) +'$'
    if sumval != 0:
      count_str += ' ('
      for c2,inherited_term in enumerate(count[m_e]):
        added= False
        if sumval == count[m_e][inherited_term]:
          count_str += 'all from $'+str(el.format_matrix_element(inherited_term))+'$)'
          break
        else:
          if count[m_e][inherited_term] != 0:
            count_str += '$'+str(count[m_e][inherited_term])+'$ from $'+str(el.format_matrix_element(inherited_term))+'$'
            added=True

        if c2 != len(count[m_e])-1:
          if added==True:
            count_str += ', '
        else:
          count_str += ')'

    if c!= len(count)-1:
      count_str += ', '
    else:
      count_str += '.'
  return count_str 

def compose_TeX_both_bases(sym,TeX_matrix_complex,TeX_matrix_real,TeX_problem,TeX_expansions_complex,TeX_expansions_real):
    TeX_str = (R'\batchmode'+'\n'
           R'\documentclass[fleqn]{article}'+'\n'
           R'\usepackage{amsmath}'+'\n'
           R'\usepackage{breqn}'+'\n'
           R'\usepackage[margin=0.7in]{geometry}' + '\n'
           R'\setlength\mathindent{0pt}'+'\n'
           R'\title{VHEGEN: A vibronic Hamiltonian expansion generator for trigonal and tetragonal polyatomic systems}'+'\n'
           R'\author{Robert A. Lang \and Riley J. Hickman \and Tao Zeng}'+'\n'
           R'\date{}'+'\n'
           R'\begin{document}'+'\n'
           R'\maketitle'+'\n'
           R'Thank you for using \texttt{VHEGEN}, the \texttt{V}-ibronic \texttt{H}-amiltonian \texttt{E}-xpansion \texttt{GEN}-erator for trigonal and tetragonal polyatomic systems. '
           R'This is a \texttt{VHEGEN} output file compiled by \texttt{pdflatex}. '
           R'If the \texttt{VHEGEN} package was used in research resulting in a publication, please reference the article in \textit{Computer Physics Communications} which describes the program ([doi here]). '
           R'Additional information regarding the matrix element expansion process, including the independent matrix element eigenvalues, their root formulas and constraints, and their transformation to the real basis (if applicable), can be found in the \texttt{log} output file. '
           R'For questions, bugs, or comments, please contact robert.lang@mail.utoronto.ca.\\\\'+'\n'
           R'\tableofcontents'+'\n'
           R'\newpage'+'\n'
           R'\section{Vibronic interaction}'+'\n'
           R''+TeX_problem+' in $'+sym[0]+'_{'+sym[1:].lower()+R'}$'+'\n'
           R'\section{Vibronic Hamiltonian operator in the complex $E$ basis}'+'\n'
           R''+TeX_matrix_complex.replace('__','')+'\n'
           R'\section{Matrix element expansions in the complex $E$ basis}'+'\n'
           R''+TeX_expansions_complex.replace('__','')+'\n'
           R'\section{Vibronic Hamiltonian operator in the real $E$ basis}'+'\n'
           R''+TeX_matrix_real.replace('__','')+'\n'
           R'\section{Matrix element expansions in the real $E$ basis}'+'\n'
           R''+TeX_expansions_real.replace('__','')+'\n'
           R'\end{document}')
    return TeX_str

def compose_TeX(sym,TeX_matrix,TeX_problem,TeX_expansions):
    TeX_str = (R'\batchmode'+'\n'
               R'\documentclass[fleqn]{article}'+'\n'
               R'\usepackage{amsmath}'+'\n'
               R'\usepackage{breqn}'+'\n'
               R'\usepackage[margin=0.7in]{geometry}' + '\n'
               R'\setlength\mathindent{0pt}'+'\n'
               R'\title{VHEGEN: A vibronic Hamiltonian expansion generator for trigonal and tetragonal polyatomic systems}'+'\n'
               R'\author{Robert A. Lang \and Riley J. Hickman \and Tao Zeng}'+'\n'
               R'\date{}'+'\n'
               R'\begin{document}'+'\n'
               R'\maketitle'+'\n'
               R'Thank you for using \texttt{VHEGEN}, the \texttt{V}-ibronic \texttt{H}-amiltonian \texttt{E}-xpansion \texttt{GEN}-erator for trigonal and tetragonal polyatomic systems. '
               R'This is a \texttt{VHEGEN} output file compiled by \texttt{pdflatex}. '
               R'If the \texttt{VHEGEN} package was used in research resulting in a publication, please reference the article in \textit{Computer Physics Communications} which describes the program. '
               R'Additional information regarding the matrix element expansion process, including the independent matrix element eigenvalues, their root formulas and constraints, and their transformation to the real basis (if applicable), can be found in the \texttt{log} output file. '
               R'For questions, bugs, or comments, please contact robert.lang@mail.utoronto.ca.\\\\'+'\n'
               R'\tableofcontents'+'\n'
               R'\newpage'+'\n'
               R'\section{Vibronic interaction}'+'\n'
               R''+TeX_problem+' in $'+sym[0]+'_{'+sym[1:].lower()+R'}$'+'\n'
               R'\section{Vibronic Hamiltonian operator}'+'\n'
               R''+TeX_matrix.replace('__','')+'\n'
               R'\section{Matrix element expansions}'+'\n'
               R''+TeX_expansions.replace('__','')+'\n'
               R'\end{document}')
    return TeX_str

def format_constraint(constraint): #{'a%ee': 'nz','broe&bieo': 'nz'}
    string = ''
    if constraint == {}:
        return ''
    for c,con in enumerate(constraint):
        cons = con.split('&')
        formatted_cons = []
        for const in cons:
            formatted_con = const
            if 'nz' in constraint[con] and con not in ['cos','sin']:
                add = 0
                if const[1] != '%':
                    formatted_con = const[:1] + '^' + const[1:]
                    add += 1
                if const[2] != '%':
                    formatted_con = formatted_con[:(2+add)] + '_' + formatted_con[(add+2):]
            formatted_cons.append(formatted_con)
        formatted_constraints = ' '.join(formatted_cons)
        formatted_constraints = formatted_constraints.replace('%','')
        string += formatted_constraints +' '+ constraint[con]
        if c != len(constraint)-1:
            string += ', '
    return string

def format_constraints(constraints):
    string = 'Constraints:\n\n'
    for e in constraints:
        if constraints[e] == {}:
            formatted_constraint = 'all nr'
        else:
            formatted_constraint = format_constraint(constraints[e])
        string += 'H_'+str(e) +' : ' + formatted_constraint +'\n'
    string +='\n'
    return string
    
def new_file(filename,extension,path):
    if path.endswith('/') != True:
        path += '/'
    file = open(path+filename+'.'+extension, "w")
    return file

def TeX_write(filename,path,TeX_str):
    TeX_file = new_file(filename,'tex',path)
    TeX_file.write(TeX_str)
    TeX_file.close()

def log_write(filename,path):
    #write log file to outputs
    global log_str
    log_file = new_file(filename,'log',path)
    log_file.write(log_str)
    log_file.close()

def convert_syntax(sympy_expr,syntax):
    if syntax.upper() not in ['LATEX','MATHEMATICA']:
        raise Exception('ConvertSyntaxError: could not recognize syntax.')
    else:
        if syntax.upper() == 'LATEX':
            return convert_to_LaTeX(sympy_expr)
        elif syntax.upper() == 'MATHEMATICA':
            return convert_to_Mathematica(sympy_expr)

def convert_to_LaTeX(sympy_expr):
    latex_expr =  sym.printing.latex(sympy_expr,fold_short_frac=True)
    latex_expr = latex_expr.replace('1.0','')
    latex_expr = latex_expr.replace('.0','')
    return latex_expr

def exec_pdflatex(filename,path):
    if not path.endswith('/'):
        path += '/'
    print('Sending output to '+path)
    if os.path.isfile(path+filename + '.pdf') == True:
        print('Overwriting '+ path+filename + '.pdf')
        os.remove(path+filename + '.pdf')
    print('Generating output PDF:\n')
    try:
        proc = subprocess.Popen(['pdflatex --interaction=batchmode -output-directory='+path, path+filename+'.tex'])
        proc.communicate()
        proc.communicate()
    except OSError:
        try:
            os.system('pdflatex --interaction=batchmode -output-directory='+path+' '+path+filename +'.tex')
            os.system('pdflatex --interaction=batchmode -output-directory='+path+' '+path+filename +'.tex') #double compile for ToC
        except OSError:
            print('Error making '+filename+'.pdf.\n')   
    '''Remove unwanted LaTeX output files .aux and .log'''
    if os.path.isfile(path+filename + '.aux') == True:
        os.remove(path+filename + '.aux')
    if os.path.isfile(path+filename + '.log') == True:
        os.remove(path+filename + '.log')
#EOF