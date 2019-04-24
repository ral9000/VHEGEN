from argparse import ArgumentParser
import modules.glbls as glo

try:
    input = raw_input
except NameError:
    pass

def configure_parser():
    parser = ArgumentParser()
    parser.add_argument("--c",dest="c",
                        help="config file path, default = config.cfg",
                        metavar="CONFIG", default="config.cfg")
    parser.add_argument("--sym",dest="sym",
                        help="point group static input",
                        metavar="SYMMETRY")
    parser.add_argument("--states",dest="states",
                        help="electronic states static input",
                        metavar="STATES")
    parser.add_argument("--modes",dest="modes",
                        help="vibrational modes static input.",
                        metavar="MODES")
    parser.add_argument("--o",dest="o",
                        help="order(s) of expansion static input.",
                        metavar="ORDERS")
    parser.add_argument("--f",dest="f",
                        help="output filename. default = 'output'",
                        metavar="FILENAME", default="output")
    return parser

def read_config(path):
    configfile = open(path,"r").readlines()
    config = {}
    for line in configfile:
        line = line.replace(u'\n','')
        if line[0] == '#':
            pass
        else:
            arg = line.split('=')
            config[arg[0]] = arg[1]
    return config

def symmetry_arg(arg):
    symmetry = arg.upper()
    if symmetry == 'LIST':
        print('\nAvailable point groups:')
        for s in glo.all_groups:
            print(s)
        raise Exception()
    elif symmetry == 'EXIT':
        exit()
    elif symmetry in glo.all_groups:
        return symmetry
    else:
        raise Exception('InputError: Symmetry '+symmetry+' not recognized.')

def states_arg(arg,sym):
    arg = arg.upper()
    arg = arg.replace('"',"''")
    if arg == 'LIST':
        print('\nIrreps of '+sym+':')
        for e in glo.irrep_dct[sym]:
            print(e)
        raise Exception()
    elif arg == 'EXIT':
        exit()
    if ',' in arg:
        states = arg.split(',')
    elif '+' in arg:
        states = arg.split('+')
    else:
        states = [arg]
    if len(states) not in (1,2):
        raise Exception('InputError: Number of states must be 1 or 2.')
    max_cond = len(states)
    met_cond = 0
    for e in states:
        if e in glo.irrep_dct[sym]:
            met_cond += 1
    if met_cond == max_cond:
        return reorder(states)
    else:
        raise Exception('InputError: State(s) '+arg+' not valid irreps in '+str(sym)+'.')

def modes_arg(arg,sym):
    arg = arg.upper()
    arg = arg.replace('"',"''")
    if arg == 'LIST':
        print('\nIrreps of '+sym+':')
        for e in glo.irrep_dct[sym]:
            print(e)
        raise Exception()
    elif arg == 'EXIT':
        exit()
    if ',' in arg:
        modes = arg.split(',')
    elif '+' in arg:
        modes = arg.split('+')
    else:
        modes = [arg]
    if len(modes) not in (1,2):
        raise Exception('InputError: Number of modes must be 1 or 2.')
    max_cond = len(modes)
    met_cond = 0
    for e in modes:
        if e in glo.irrep_dct[sym]:
            met_cond += 1
    if met_cond == max_cond:
        return reorder(modes)
    else:
        raise Exception('InputError: Modes(s) '+modes.lower()+' not valid irreps in '+symmetry+'.')

def orders_arg(arg):
    if arg.upper() == 'EXIT':
        exit()
    orders = arg.split(',')
    for i,order in enumerate(orders):
        try:
            orders[i] = int(order)
        except ValueError:
            raise Exception('InputError: Orders must be non-negative integers.')
    if len(orders) > 2:
        raise Exception('InputError: Invalid range of orders.')
    elif len(orders) == 2:
        if orders[1] <= orders[0]:
            raise Exception('InputError: Invalid range of orders.')
        if (orders[0]) < 0 or (orders[1] < 0):
            raise Exception('DynamicInputError: Orders must be non-negative integers.')
    elif len(orders) == 1:
        if orders[0] < 0:
            raise Exception('InputError: Orders must be non-negative integers.')
    if len(orders) == 2:
        orders = range(int(orders[0]),int(orders[1])+1)
    return orders
    
def filename_arg(arg):
    for char in (r"!@#$%^&*(){}[],."):
        if char in arg:
            print(char)
            raise Exception('InputError: Filename may only contain letters, numbers, "_", and "-".')
    if len(arg) == 0:
        return 'output'
    else:
        return arg

def format(inp_list):
    formatted_str = '('
    for i,n in enumerate(inp_list):
        if i==1:
            formatted_str+='+'
        formatted_str+=n
    formatted_str+=')'
    return formatted_str

def format_problem(states,modes):
    states_str = format(states)
    modes_str = format(modes)
    return states_str+'x'+modes_str.lower()

def reorder(inp_list):
    #re-orders states/modes by priority eg) (A+E)x(e+a) -> (E+A)x(e+a)
    if len(inp_list) > 1:
        if glo.irrep_priority.index(inp_list[0][0]) > glo.irrep_priority.index(inp_list[1][0]):
            inp_list = list(reversed(inp_list))
    return inp_list
    
def return_problem(symmetry,states,modes,orders,filename):
    return  ('\n---------------------------------------'
             '\n  VHEGEN instance parameters:'
             '\n  Symmetry: ' + symmetry+
             '\n  Problem: ' + format_problem(states,modes)+
             '\n  Expansion orders: ' + str(orders)+
             '\n  Output filename: ' + filename+
             '\n---------------------------------------')

def dynamic_input():
    #symmetry
    dyn_inp = input('Enter symmetry: ')
    while True:
        try:
            sym_inp = symmetry_arg(dyn_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Re-enter symmetry: ')
            continue
        else:
            symmetry = sym_inp
            print(symmetry +' symmetry accepted.')
            break
    #states
    dyn_inp = input('Enter electronic state(s): ')
    while True:
        try:
            states_inp = states_arg(dyn_inp,symmetry)
            states_inp = reorder(states_inp)
            if symmetry in glo.trigonal_groups:
                states_inp = states_inp #limit_trigonal(states_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Retry state(s): ')
            continue
        else:
            states = states_inp
            print('Electronic states '+ format(states)+' accepted.')
            break
    #modes
    dyn_inp = input('Enter vibrational mode(s): ')
    while True:
        try:
            modes_inp = states_arg(dyn_inp,symmetry)
            modes_inp = reorder(modes_inp)
            if symmetry in glo.trigonal_groups:
                modes_inp = modes_inp #limit_trigonal(modes_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Retry mode(s): ')
            continue
        else:
            modes = modes_inp
            print('Vibrational modes '+format(modes).lower()+' accepted.')
            break
    #orders
    dyn_inp = input('Order(s) of expansion:')
    while True:
        try:
            orders_inp = orders_arg(dyn_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Retry orders: ')
            continue
        else:
            orders = orders_inp
            print('Orders of expansion '+(str(list(orders)))+' accepted.')
            break
    #output filename
    dyn_inp = input('Enter filename: ')
    while True:
        try:
            filename_inp = filename_arg(dyn_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Retry filename: ')
            continue
        else:
            filename = filename_inp
            print('Filename: "'+filename+'".')
            break

    print(return_problem(symmetry,states,modes,orders,filename))

    while True:
        yn_prompt = input('Continue? [Y/N]: ')
        if yn_prompt.upper() == 'Y':
            print('')
            return symmetry, states, modes, orders, filename
        elif yn_prompt.upper() == 'N' or yn_prompt.upper() == 'EXIT':
            exit()

def read_input():
    argparser = configure_parser()
    options = argparser.parse_args()
    config = read_config('config.cfg')
    if config[u'input'] == u'static':
        print("Reading static inputs.")
        #static input
        try:
            symmetry = symmetry_arg(options.sym)
        except Exception as e:
            print('StaticInputError: Invalid symmetry argument.')
            exit()
        try:
            states = reorder(states_arg(options.states,symmetry))
            if symmetry in glo.trigonal_groups:
                states = states #limit_trigonal(states)
        except Exception as e:
            print('StaticInputError: Invalid states argument.')
            exit()
        try:
            modes = reorder(modes_arg(options.modes,symmetry))
            if symmetry in glo.trigonal_groups:
                modes = modes #limit_trigonal(modes)
        except Exception as e:
            print(e)
            print('StaticInputError: Invalid modes argument.')
            exit()
        try:
            orders = orders_arg(options.o)
        except Exception as e:
            print('StaticInputError: Invalid orders argument.')
            exit()
        try:
            filename = filename_arg(options.f)
        except Exception as e:
            print('StaticInputError: Invalid filename argument.')
            exit()
    elif config[u'input'] == u'dynamic':
        print("Entering dynamic input.")
        #enter dynamic input
        symmetry,states,modes,orders,filename = dynamic_input()
    problem_dct = {'sym':symmetry,
                   'states':states,
                    'modes':modes,
                    'o': orders,
                    'f': filename}
    else:
        print("Error: 'input' value in config.cfg not found or recognized.")
        exit()
    return config,problem_dct

def prepare_input(sym,states,modes,orders,filename='output'):
    sym = symmetry_arg(sym)
    states = states_arg(states,sym)
    modes = modes_arg(modes,sym)
    orders = orders_arg(orders)
    filename = filename_arg(filename)
    return {'sym':sym,
            'states':states,
            'modes':modes,
            'o': orders,
            'f':filename}
#EOF