#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Run ST-on-the-fly"""

#------------------------------------------------------------------------------
# PEP8 imports : 
# on separate lines
# at the top of the file


# 1. standard library imports
#---
import argparse
from __future__ import print_function
from __future__ import division 
import math
import os
import pdb
import random
import shlex
import shutil
import subprocess
import sys

# 2. related third party imports
#---
import numpy as np
# scipy is a library (package) that contains modules 
# and to import a specific module from the scipy library, 
# it is needed to specify it and import the module itself.
# So do not : 
# >>> import scipy
# >>> scipy.constants.<etc>
# but : 
from scipy import constants


# 3. local application/library specific imports
#---
import logger



def get_arguments_values(): 
    """ Use argparse module to get arguments values.
    Return
    ------
    args : Namespace object
        Namespace object populated with arguments converted from strings to objects 
        and assigned as attributes of the namespace object.
        They store arguments values (str or int or whatever, according to the specifications in
        the code)
    """
    help_str = "ST-on-the-fly experiment"

    parser = argparse.ArgumentParser(description=help_str,add_help=True)

    # argparse automatically convert any internal - characters to _ characters 
    # to make sure the string is a valid attribute name.
    # So in this case, use dest is not necessary.

    # About T_RANGE : 
    parser.add_argument("--Tmin", required=True,
        help="At which temperature (in Kelvin) begins the experiment", type=float)
    parser.add_argument("--Tmax", required=True,
        help="The maximal temperature (in Kelvin) during the experiment", type=float)
    parser.add_argument("--Tstep", required=True,
        help="The step of the range of temperature ", type=int)


    # About struct : 

    parser.add_argument("--gro-filename", required=True,
        help="gro / pdb file to use for the ST experiment", type=str)
    parser.add_argument("--top-filename", required=True,
        help="topology file to use for the ST experiment ", type=str)

    # About ST : 
    parser.add_argument("--nb-md",required=True,
        help="The number of molecular dynamics during the experiment : t_one-md * num-md = t_ST", type=int)

    parser.add_argument("--st-mdp-template-filename", required=True,
        help="mdp file to use for the ST experiment", type=str)

    parser.add_argument("--st-outname", required=True,
        help="template name for output of the ST experiment ", type=str)


    # About minimisation : 
    minimisation_parser = parser.add_mutually_exclusive_group(required=False)
    minimisation_parser.add_argument('--minimisation', 
        action='store_true', 
        help="Run a minimisation before ST experiment")
    minimisation_parser.add_argument('--no-minimisation', 
        dest='minimisation', action='store_false', 
        help="Do run a minimisation before ST experiment")
    parser.set_defaults(minimisation=False)

    parser.add_argument("--minimisation-mdp-filename", 
        help="mdp file to use for minimisation", type=str, required=False)
    parser.add_argument("--minimisation-outname", 
        help="template name for output of minimisation ", type=str, required=False)


    # Other :
    parser.add_argument("--out-path", 
        help = "Where the outputed results files should be store", type= str, default='./')

    parser.add_argument("--maxwarn", 
        help="The max number of warnigs allowed when running MD", type=int, default = '0')

    parser.add_argument("--clean-all", 
        help="Clean gromacs outputs unecessary to plot the figure in N'Guyen 2013 /!\ Be carefull", 
        action='store_true')

    parser.add_argument('-v', '--verbose', dest = 'verbosity',
        help="Turn on detailed info log", action='count',  default=0)

    return parser.parse_args()


def print_use_help_message() : 
    logger.__logger.error("Use --help for help.\n")


def assert_strictly_positive(value, name_var): 
    try:
        assert (value > 0 )
    except AssertionError, e:
        logger.__logger.error(name_var + ' must be strictly positive.')
        raise e 

def assert_strictly_inferior(value1, value2, name_var1, name_var2): 
    try:
        assert (value1 < value2 )
    except AssertionError, e:
        logger.__logger.error(name_var1 + ' must be strictly inferior to ' + name_var2+'.')
        raise e 


def check_arguments_integrity(args): 
    """ Check integrity of arguments pass through the command line. 
    Parameter
    ---------
    args : namespace object
        Its attributes are arguments names 
        and contain arguments values. 
    Side effect
    -----------
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
        assert_strictly_inferior (args.Tmin + args.Tstep,  args.Tmax, 'Tmin + Tstep', 'Tmax') 
    except AssertionError:
        print_use_help_message()
        sys.exit(-1)


    try:
        assert_strictly_positive(args.nb_md, 'nb-md')
    except AssertionError:
        print_use_help_message()
        sys.exit(-1)


    files = [args.gro_filename, args.top_filename, args.st_mdp_template_filename]
    if args.minimisation : 
        files.append (args.minimisation_mdp_filename)

    for file in files : 
        try:
            with open(file):
                pass
        except IOError as e:
            logger.__logger.error ("Unable to open file "+ file +". (Does not exist or no read permissions)" ) 


def get_integrous_arguments_values(): 
    """ Get arguments passed through command line and check their integrity.
    Return
    ------
    args : namespace object
        Its attributes are arguments names 
        and contain arguments values (str or int or whatever according to code specifications). 
    """
    args = get_arguments_values()
    check_arguments_integrity(args)
    return args


class ListWithoutNegIdx(list):

    def __getitem__(self, key):
        if isinstance(key, int):
            if key < 0:
                raise IndexError("negative index not allowed")
            else :  
                return super(ListWithoutNegIdx, self).__getitem__(key)
        else : 
            raise TypeError ("ListWithoutNegIdx indices must be integers, not "+ str(type (key) ) )  




class MolecularDynamics(object):
    """docstring for MolecularDynamics"""
    @logger.log_decorator
    def __init__(self, mdp_filename, gro_filename, top_filename, out_path, out_name, maxwarn, **kwargs):
        print ('je suis dans le constructor de MolecularDynamics')
        super(MolecularDynamics, self).__init__(**kwargs)
        self._mdp_filename = mdp_filename
        self._out_path = out_path
        self._out_name = out_name
        self._gro_filename = self.out_path + self.out_name + '.gro'
        shutil.copyfile(gro_filename, self._gro_filename)
        self._top_filename = top_filename
        self._maxwarn = maxwarn

    @property 
    def mdp_filename(self):
        return self._mdp_filename
    @property 
    def gro_filename(self):
        return self._gro_filename
    @property 
    def top_filename(self):
        return self._top_filename
    @property 
    def out_path(self):
        return self._out_path
    @property 
    def out_name(self):
        return self._out_name
    @property 
    def maxwarn(self):
        return self._maxwarn

    @logger.log_decorator
    def run(self):
        print ('MD.run()=========')

    @logger.log_decorator
    def gmx_grompp(self) : 
        cmd = "gmx grompp -f {0} -c {1} -p {2} -o {3}.tpr -po {3}.mdp -maxwarn {4}".format(   \
            self.mdp_filename, \
            self.gro_filename, \
            self.top_filename, \
            self.out_path + self.out_name , \
            self.maxwarn \
        )
        print (cmd)
        p = subprocess.Popen(shlex.split(cmd))
        p.wait()

    @logger.log_decorator
    def gmx_mdrun(self):
        cmd = "gmx mdrun -v -s {0}.tpr -o {0}.trr -e {0}.edr -g {0}.log -c {0}.gro".format(self.out_path + self.out_name)
        print (cmd)
        p = subprocess.Popen(shlex.split(cmd))
        p.wait()

    @logger.log_decorator
    def run(self): 
        self.gmx_grompp()
        self.gmx_mdrun()

class Simulation(object):
    """docstring for Simulation"""
    @logger.log_decorator
    def __init__(self, T_current, **kwargs):
        print ('je suis dans le constructor de Simulation')

        super(Simulation,self).__init__(**kwargs)
        print ('je suis revenu de super')
        self._T_current = T_current

    @property
    def T_current(self):
        return self._T_current
    


class MonteCarlo(Simulation):
    """docstring for MonteCarlo"""
    
    @logger.log_decorator
    def run(self):
        return compute_E_average()
    @logger.log_decorator
    def compute_E_average(self):
        return 0

class MolecularDynamicsProduction(Simulation,MolecularDynamics):
    """docstring for MolecularDynamicProduction"""
    @logger.log_decorator
    def __init__(self, **kwargs):
        print ('je suis dans le constructor de Prod') 
        super(MolecularDynamicsProduction, self).__init__(**kwargs)
        self._mdp_template = kwargs['mdp_filename']
        self._mdp_filename = '{0}_{1}.mdp'.format(self._mdp_template , self.T_current)


    @logger.log_decorator
    def run(self):
        super(MolecularDynamicsProduction, self).run() # call MolecularDynamics.run()
        return self.compute_E_average()

    @logger.log_decorator
    def gmx_energy(self, arg = 'Potential') : 

        p1 = subprocess.Popen( shlex.split("echo {0}".format(arg)), stdout=subprocess.PIPE ) 

        cmd = "gmx energy -f {0}.edr -o {0}_Potential.xvg ".format(
                    self.out_path + self.out_name
                )
               
        p2 = subprocess.Popen(shlex.split(cmd), stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        p2_com = p2.communicate() 
        p2_out = p2_com[0].split('\n') 

        return p2_out

    @logger.log_decorator
    def compute_E_average(self):
        output = self.gmx_energy('Potential')
        for line in output : 
            splitted_line = line.split()
            print (splitted_line)
            if splitted_line and splitted_line[0] == 'Potential' : 
                E_average = splitted_line[1]
                return int(E_average)
        log.error('E_average could not be found') 
        return 0

    @Simulation.T_current.setter
    def T_current(self, T_new):
        pdb.set_trace()
        print ("setter de Tcurrent ============")
        self.update_velocities(T_new)

        self._T_current = T_new 
        # not self.T_current = T_new
        # otherwise it calls the setter that calls the setter, that c... 
        # do not use setter in setter ! 
        
        self._mdp_filename = self._mdp_template % self.T_current

        
        

    def update_velocities(self, T_new):

        #Can't change a file, so create a new temp one...
        nb_of_atom = 0 

        with open(self.gro_filename, 'r') as infile, \
            open(self.gro_filename+'.tmp', 'w') as outfile    :
            for line_idx, line in enumerate (infile) : 

                # gro file line : 
                #"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" + \n
                # |___________44_____________|              |
                # |____________________68___________________|

                # nb of atom on seconde line -> 2-1
                if line_idx == 1 : 
                    nb_of_atom = int(line )

                # from third line to x = 3 + nb_of_atom th line : 
                if line_idx in xrange(2, 2 + nb_of_atom + 1) : 
                    #Remember [a, b] = xrange (a, b + 1 )

                    if len(line) < 69: # 68 + \n 
                        log.error('There is no velocities in the {0} file'.format(infile))
                        exit(-1)

                    velocities = [ float(line[44:52]), float(line[52:60]), float(line[60 :68 ]) ]
                    new_velocities = list()
                    for velocity in velocities : 
                        new_velocity = velocity * math.sqrt(T_new / self.T_current)
                        new_velocities.append(new_velocity)
                    
                    line = "{0}{1:8.4f}{2:8.4f}{3:8.4f}\n".format(line [:44], *new_velocities)
                    print (line)

                outfile.write(line)
        print ("làààààààààààààààààààààààààààààààààààààààà")   
        # ... and rename it (the old one is erase)
        os.rename(self.gro_filename+'.tmp', self.gro_filename)



    
@logger.log_decorator
def create_simulation(simu_type, **kwargs): 
    if simu_type == 'md' : 
        return MolecularDynamicsProduction(**kwargs)
    elif simu_type == 'mc':
        return MonteCarlo(**kwargs)
    else : 
        logger.__logger.error('Wrong choice of simulation')
        return None



class NoECurrent(Exception):
    def __init__(self):
        super(NoECurrent, self).__init__()

class Temperature(object):
    """docstring for Temperature"""

    k_Boltzmann = constants.value(u'Boltzmann constant') # ??? k

    def __init__(self, value):
        super(Temperature, self).__init__()
        self._VALUE = value
        self._number_of_passes = 0 
        self._f = 0
        self._E = None
        self._BETA  = self.compute_beta()


    @logger.log_decorator
    def compute_beta(self):
        return 1/Temperature.k_Boltzmann * self._VALUE
        # __future__ division -> floting point division


    @logger.log_decorator
    def update_f(self, Tprev) : 
        try:
            self.compute_f(Tprev)
        except NoECurrent:
            self.estimate_f(Tprev)


    @logger.log_decorator
    def estimate_f(self, Tprev):
        # (beta_next - beta_previous) E_previous / 2
        self._f = (self._BETA - Tprev._BETA ) * Tprev._E / 2


    @logger.log_decorator
    def compute_f(self, Tprev):
        # f_prev + (beta_curr - beta_prev) (E_curr + E_prev) / 2 
        try:
            self._f =  Tprev._f + (self._BETA - Tprev._BETA ) * ( self._E + Tprev._E )  / 2
        except Exception:
            raise NoECurrent


    @logger.log_decorator
    def update_E(self, E_new):
        self._number_of_passes += 1 
        try:
            self._E =  self._E  + ( (E_new - self._E ) / self._number_of_passes)
        except TypeError: 
            #self._E = None becase no updated yet
            #First time init : 
            self._E =  ( E_new  / self._number_of_passes)
        
        


class SimulatedTempering(object):
    """docstring for ST"""

    @logger.log_decorator
    def __init__(self, num_simu, Tmin, Tmax, Tstep, simu_type='md', st_mdp_template_filename = None, **kwargs):
        
        super(SimulatedTempering,self).__init__()
        self._NUM_SIMU = num_simu
        self._T_RANGE=ListWithoutNegIdx() 
        self._ST_MDP_TEMPLATE_FILENAME= st_mdp_template_filename

        for T in np.arange(Tmin,Tmax+1,Tstep):
            if simu_type == 'md' : 
                self.create_mdp(T)
            self._T_RANGE.append(Temperature(T))
        #range (a, b) = [a, b[
        #range (a, b+1) = [a, b+1[ = [a, b]
        kwargs['T_current'] = self._T_RANGE[0]._VALUE
        kwargs['mdp_filename'] = self._ST_MDP_TEMPLATE_FILENAME
        self._SIMULATION=create_simulation(simu_type, **kwargs ) #pattern strategy



    @logger.log_decorator
    def create_mdp(self, T) : 
        
        mdp_filename = '{0}_{1}.mdp'.format(self._ST_MDP_TEMPLATE_FILENAME, T)

        with open(self._ST_MDP_TEMPLATE_FILENAME, 'r') as infile, \
            open(mdp_filename, 'w') as outfile    :
            for line in infile : 
                if line.startswith('ref_t') : 
                    line = "ref_t = {0}\n".format(T)

                outfile.write(line)


    @property
    def T_current(self):
        return self._T_RANGE[self.T_current_idx]
    
    @property    
    def T_current_idx(self):
        return self.get_T_idx(self._SIMULATION.T_current)

    def get_T_idx(self, T_wanted) : 
        return [T._VALUE for T in self._T_RANGE].index(T_wanted)
        
    @logger.log_decorator
    def update_f_current(self):
        try : 
            T_previous = self._T_RANGE[self.T_current_idx-1]
            self.T_current.update_f(T_previous)
        except IndexError : #no previous T because Tcurrent = Tmin 
            self.T_current._f = 0 #f_Tmin is always equal to 0.

        # Remember : 
        # If an exeption occurs during execution of the try clause,
        # the rest of the clause is skipped.    

    @logger.log_decorator
    def update_f_next(self):
        try : 
            T_next = self._T_RANGE[self.T_current_idx + 1 ]
            T_next.update_f( self.T_current)
        except IndexError : #no next T because Tcurrent = Tmax 
            pass

        # Remember : 
        # If an exeption occurs during execution of the try clause,
        # the rest of the clause is skipped.    



    # @staticmethod returns descriptor objects, not functions. 
    # problem : most decorators are not designed to accept descriptors.
    # solution : @staticmethod must be the top-most decorator 
    # for another decorator to decorate it.
    @staticmethod
    @logger.log_decorator
    def toss_coin():
        return random.choice ( [-1, 1] )


    @logger.log_decorator
    def choose_T_attempt(self):
        try:
            T_attempt = self._T_RANGE[self.T_current_idx + SimulatedTempering.toss_coin() ]
        except IndexError:
            T_attempt = None
        return T_attempt

    @logger.log_decorator
    def compute_metropolis_criterion(self, T_attempt) : 

        try:
            insider = - (                                                        \
                ( T_attempt._BETA- self.T_current._BETA )  * self.T_current._E   \
                - ( T_attempt._f - self.T_current._f )                           \
            )           
            res = math.exp ( insider  ) 
        except OverflowError: 
            if insider < 0 : 
                res = 0
            elif insider > 0 : 
                res = float('inf')
            else : 
                print ('case not handled !! ')
                pdb.set_trace()
        
        mc = min (1, res)
        return mc
        # mc = min (1 , exp (-   [  (beta_attempt - beta_current) * E_current_average  -   (f_attempt_estimate - f_current_compute)  ]   ) )

    @logger.log_decorator
    def attempt_OK(self, T_attempt):
        if not T_attempt : return False
        mc =  self.compute_metropolis_criterion(T_attempt) 
        if mc == 1 : 
            return True
        else : 
            return False

    @logger.log_decorator
    def run(self):
        for step_idx in xrange(self._NUM_SIMU) : 
            pdb.set_trace()
            E_current_average = self._SIMULATION.run()
            pdb.set_trace()
            self.T_current.update_E(E_current_average) 

            self.update_f_current()
            self.update_f_next ()


            T_attempt = self.choose_T_attempt()
        
            if self.attempt_OK(T_attempt) : 
                self._SIMULATION.T_current = T_attempt._VALUE #if MD : change velocity --> Overriding
            else : 
                self._SIMULATION.T_current = self.T_current._VALUE #if MD : change velocity --> Overriding



if __name__ == "__main__":
    """
    python ST-on-the-fly.py --Tmin 1 --Tmax 5 --Tstep 1 \
    --gro-filename ../data/ala10_md000.pdb \
    --top-filename ../data/ala10.top \
    --nb-md 3 \
    --st-mdp-template-filename ../data/md1.mdp \
    --st-outname outst \
    --minimisation \
    --minimisation-mdp-filename ../data/mini2.mdp  \
    --minimisation-outname miniout \
    --maxwarn 1 \
    --out-path ./ \
    -vvv
    """
    print ('go')
    logger.set_functional_logger()

    logger.__logger.info('logger OK')

    logger.__logger.info('get args')
    args = get_integrous_arguments_values()
    


    verbosity = [ 0, 30  ,20 ,10]
    try:
        level = verbosity[args.verbosity] 
    except IndexError:
        level = 10 

    logger.__logger.setLevel(level)

    global log
    log = logger.__logger
    if args.minimisation : 

        logger.__logger.info('run minimi')
        MD_minimi = MolecularDynamics(        
            mdp_filename = args.minimisation_mdp_filename, 
            gro_filename = args.gro_filename, 
            top_filename = args.top_filename, 
            out_path = args.out_path, 
            out_name = args.minimisation_outname, 
            maxwarn = args.maxwarn)
        MD_minimi.run()
        logger.__logger.info('minimi OK')

        st_gro_filename = args.out_path + args.minimisation_outname + '.gro'
    else : 
        st_gro_filename = args.gro_filename
    logger.__logger.info('new ST')
    ST = SimulatedTempering(
        args.nb_md, 
        args.Tmin, args.Tmax, args.Tstep, 
        'md',  
        st_mdp_template_filename = args.st_mdp_template_filename, 
        gro_filename = st_gro_filename, 
        top_filename = args.top_filename, 
        out_path = args.out_path, 
        out_name = args.st_outname, 
        maxwarn = args.maxwarn
    )
    
    logger.__logger.info('ST.run')
    ST.run()

    logger.__logger.info('the end')