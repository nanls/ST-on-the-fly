#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Run ST-on-the-fly

Example : 
---------

# With previous minimisation : 
python ST-on-the-fly.py \
--Tmin 1 --Tmax 5 --Tnum 1 \
--gro-filename ../data/ala10_md000.gro \
--top-filename ../data/ala10.top \
--nb-md 3   --st-mdp-template-filename ../data/md1.mdp  --st-outname outst  \
--minimisation \
--minimisation-mdp-filename ../data/mini2.mdp  --minimisation-outname miniout \ 
--maxwarn 1     --out-path ./     -vvv

#Without previous minimisation : 
python ST-on-the-fly.py \
--Tmin 1 --Tmax 5 --Tnum 1  \
--gro-filename ../data/ala10_md000.gro \
--top-filename ../data/ala10.top  \
--nb-md 3  --st-mdp-template-filename ../data/md1.mdp  --st-outname outst \
--maxwarn 1     --out-path ./     -v

"""

#------------------------------------------------------------------------------
# PEP8 imports : 
# on separate lines
# at the top of the file

# 1. standard library imports
#---

# from __future__ imports must occur at the beginning of the file.
# if not : SyntaxError 
from __future__ import print_function
from __future__ import division 

import argparse
import errno
import pdb
import os

# 2. other imports
#--- 

# 3. local imports
#---
import logger
logger.set_functional_logger()
global log
log = logger.__logger
log.info('logger OK')

from md import MolecularDynamics
from st import SimulatedTempering


#------------------------------------------------------------------------------
# Arguments getter functions : 

def get_arguments_values(): 
    """ Use argparse module to get arguments values.

    Return
    ------
    args : Namespace object
        Namespace object populated with arguments converted from strings to 
        objects (str or int, or whatever, according to the spec in the code)
        and assigned as attributes of the namespace object.
    """
    help_str = "ST-on-the-fly experiment"

    parser = argparse.ArgumentParser(description=help_str,add_help=True)

    # argparse automatically convert any internal - characters to _ characters 
    # to make sure the string is a valid attribute name.
    # So in this case, use dest is not necessary.

    # About T_range : 
    parser.add_argument("--Tmin", required=True, type=float,
        help="At which temperature (in Kelvin) begins the experiment")
    parser.add_argument("--Tmax", required=True, type=float,
        help="The maximal temperature (in Kelvin) during the experiment")
    parser.add_argument("--Tnum", required=True, type=int,
        help="Number of temperature in the given range [Tmin, Tmax]")

    # About struct : 
    parser.add_argument("--gro-filename", required=True, type=str,
        help="gro / pdb file to use for the ST experiment")
    parser.add_argument("--top-filename", required=True, type=str,
        help="topology file to use for the ST experiment " )

    # About ST : 
    parser.add_argument("--nb-md",required=True, type=int,
        help=("The number of molecular dynamics during the experiment :"
            "t_one-md * num-md = t_ST")
    )
    parser.add_argument("--st-mdp-template-filename", required=True, type=str,
        help="mdp file to use for the ST experiment")
    parser.add_argument("--st-outname", required=True, type=str, 
        help="template name for output of the ST experiment ")

    # About minimisation : 
    minimisation_parser = parser.add_mutually_exclusive_group(required=False)
    minimisation_parser.add_argument('--minimisation', 
        action='store_true', 
        help="Run a minimisation before ST experiment")
    minimisation_parser.add_argument('--no-minimisation', 
        dest='minimisation', action='store_false', 
        help="Do run a minimisation before ST experiment")
    parser.set_defaults(minimisation=False)
    parser.add_argument("--minimisation-mdp-filename", required=False,type=str,
        help="mdp file to use for minimisation")
    parser.add_argument("--minimisation-outname", required=False, type=str,
        help="template name for output of minimisation"  )
    parser.add_argument("--gene_veloc_outname", required=False, type=str,
        help="template name for output of velocities generation"  )

    # Other :
    parser.add_argument("--out-path",default='./',  type= str, 
        help = "Where the outputed results files should be store")
    parser.add_argument("--maxwarn", default = '0', type=int, 
        help="The max number of warnigs allowed when running MD" )
    parser.add_argument("--clean-all", 
        help=("Clean gromacs outputs unecessary to plot the figures"
        "in N'Guyen 2013 /!\ Be carefull"), 
        action='store_true'
    )
    parser.add_argument('-v', '--verbose', dest = 'verbosity', default=0,
        help="Turn on detailed info log", action='count')

    return parser.parse_args()


def print_use_help_message() : 
    """ Print use help message
    """
    log.error("Use --help for help.\n")


def assert_strictly_positive(value, name_var): 
    """ Assert a value is strictly positive. If not, log an error. 

    Arguments : 
    ------------
    value : number
        The value of which you want to test the strict positivity 
    name_var : string
        The name of the value

    Side effect : 
    -------------
    If the assertion fails, it logs an error using logger module:
    <name_var> must be strictly positive

    Example : 
    ---------

    >>> assert_strictly_positive (1, 'my_var_name') 
    >>> assert_strictly_positive (-1, 'my_var_name') 
    2016-09-18 15:14:17,262 :: ERROR :: my_var_name must be strictly positive.
    """
    try:
        assert (value > 0 )
    except AssertionError, e:
        log.error(name_var + ' must be strictly positive.')
        raise e 

def assert_strictly_inferior(value1, value2, name_var1, name_var2): 
    """Assert a value is strictly inferior to another one. If not, log an error. 

    Arguments : 
    ----------
    value1 : number
        The value of which you want to test the strict inferiority compared
        to value2
    value2 : number 
        The value of which you want to test the strict superiority compared 
        to value1
    name_var1 : string
        The name of value1
    name_var2 : string 
        The name of value2

    Side effect : 
    -------------
    If the assertion fails, it logs an error using logger module:
    <name_var1> must be strictly inferior to <name_var2>

    """
    try:
        assert (value1 < value2 )
    except AssertionError, e:
        log.error(name_var1 + ' must be strictly inferior to ' + name_var2+'.')
        raise e 


def check_arguments_integrity(args): 
    """ Check integrity of arguments pass through the command line. 

    Arguments :
    -----------
    args : namespace object
        Its attributes are arguments names and contain arguments values. 
    Side effect : 
    -------------
    If one arg is not integrous, exit the script printing why. 
    """
    try:
        assert_strictly_positive(args.Tmin, 'Tmin')
    except AssertionError:
        print_use_help_message()
        sys.exit(-1)
    
    try:
        assert_strictly_positive(args.Tmax, 'Tmax')
    except AssertionError:
        print_use_help_message()
        sys.exit(-1)

    try:
        assert_strictly_inferior(args.Tmin, args.Tmax, 'Tmin', 'Tmax')
    except AssertionError:
        print_use_help_message()
        sys.exit(-1)

    try:
        assert_strictly_inferior (args.Tmin + args.Tnum,  args.Tmax, 
            'Tmin + Tnum', 'Tmax') 
    except AssertionError:
        print_use_help_message()
        sys.exit(-1)


    try:
        assert_strictly_positive(args.nb_md, 'nb-md')
    except AssertionError:
        print_use_help_message()
        sys.exit(-1)


    files=[args.gro_filename, args.top_filename, args.st_mdp_template_filename]
    if args.minimisation : 
        files.append (args.minimisation_mdp_filename)

    for file in files : 
        try:
            with open(file):
                pass
        except IOError as e:
            log.error ("Unable to open file "+ file +\
                ". (Does not exist or no read permissions)" ) 


def get_integrous_arguments_values(): 
    """ Get arguments passed through command line and check their integrity.

    Return : 
    --------
    args : namespace object
        Its attributes are arguments names and contain arguments values 
        (str or int or whatever according to code specifications). 
    """
    args = get_arguments_values()
    check_arguments_integrity(args)
    return args

#------------------------------------------------------------------------------

def make_sure_path_exists(path):
    """Make sure a dir exists 

    Argument : 
    ----------
    path : string 
        The path you wan't to assert the existance 

    Raise : 
    --------
    OSError if problem when creating the dir
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise



################################################################################
if __name__ == "__main__":

    print ('go')

    #-----------------
    # get arguments : 
    log.info('get args')
    args = get_integrous_arguments_values()
    make_sure_path_exists(args.out_path)

    #-----------------
    # set verbosity : 
    level = logger.get_level(args.verbosity)
    log.setLevel(level)

    #-----------------
    # run minimi + velocities generation if needed : 
    if args.minimisation : 

        log.info('run minimi')
        MD_minimi = MolecularDynamics(        
            mdp_filename = args.minimisation_mdp_filename, 
            gro_filename = args.gro_filename, 
            top_filename = args.top_filename, 
            out_path = args.out_path, 
            out_name = args.minimisation_outname, 
            maxwarn = args.maxwarn)
        MD_minimi.run()
        log.info('minimi OK')

        st_gro_filename = args.out_path + args.minimisation_outname + '.gro'

        log.info('gene velocities')
        MD_gene_velo = MolecularDynamics(        
            mdp_filename = "../data/gene_velocities.mdp", 
            gro_filename = st_gro_filename, 
            top_filename = args.top_filename, 
            out_path = args.out_path, 
            out_name = args.gene_veloc_outname, 
            maxwarn = args.maxwarn)
        MD_gene_velo.run()
        log.info('gene vel OK')

        st_gro_filename = args.out_path + args.gene_veloc_outname + '.gro'


    else : 
        st_gro_filename = args.gro_filename

    #------------------
    # ST experiment : 
    log.info('new ST')
    ST = SimulatedTempering(
        args.nb_md, 
        args.Tmin, args.Tmax, args.Tnum, 
        'md',  
        st_mdp_template_filename = args.st_mdp_template_filename, 
        gro_filename = st_gro_filename, 
        top_filename = args.top_filename, 
        out_path = args.out_path, 
        out_name = args.st_outname, 
        maxwarn = args.maxwarn
    )
    
    log.info('ST.run')
    ST.run()

    log.info('the end')