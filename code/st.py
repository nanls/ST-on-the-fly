#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


import numpy as np
# scipy is a library (package) that contains modules 
# and to import a specific module from the scipy library, 
# it is needed to specify it and import the module itself.
# So do not : 
# >>> import scipy
# >>> scipy.constants.<etc>
# but : 
from scipy import constants


import logger
from md import MolecularDynamics


logger.set_functional_logger()
global log
log = logger.__logger


@logger.log_decorator
def create_simulation(simu_type, **kwargs): 
    if simu_type == 'md' : 
        return MolecularDynamicsProduction(**kwargs)
    elif simu_type == 'mc':
        return MonteCarlo(**kwargs)
    else : 
        logger.__logger.error('Wrong choice of simulation')
        return None

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
                if line_idx in xrange(2, 2 + nb_of_atom ) : 
                	# this range do not take the line defining the box

                    if len(line) < 69: # 68 + \n 
                        log.error('There is no velocities in the {0} file'.format(infile))
                        exit(-1)

                    velocities = [ float(line[44:52]), float(line[52:60]), float(line[60 :68 ]) ]
                    new_velocities = list()
                    for velocity in velocities : 
                        new_velocity = velocity * math.sqrt(T_new / self.T_current)
                        new_velocities.append(new_velocity)
                    
                    line = "{0}{1:8.4f}{2:8.4f}{3:8.4f}\n".format(line [:44], *new_velocities)
                    

                outfile.write(line)
        # ... and rename it (the old one is erase)
        os.rename(self.gro_filename+'.tmp', self.gro_filename)






class SimulatedTempering(object):
    """SimulatedTempering class
    

    """

    class TRange(list):
        """A list without negative indiciation
        """

        class NegativeIndexError(IndexError):
            """docstring for Exception"""
            def __init__(self,*args,**kwargs):
                super(Exception, self).__init__(*args,**kwargs)


        def __getitem__(self, key):
            """Override list getitem to disable negative indices
            """
            if isinstance(key, int):
                if key < 0:
                    
                    raise SimulatedTempering.TRange.NegativeIndexError ("neg index not allowed")
                    
                else :  
                    return super(SimulatedTempering.TRange, self).__getitem__(key)
            else : 
                raise TypeError ("TRange indices must be integers, not "+ str(type (key) ) )  




    class Temperature(object):
        """docstring for Temperature"""

        class NoECurrent(Exception):
            def __init__(self):
                super(SimulatedTempering.Temperature.NoECurrent, self).__init__()


        k_Boltzmann = constants.value(u'Boltzmann constant') # ??? k

        def __init__(self, value):
            super(SimulatedTempering.Temperature, self).__init__()
            self._VALUE = value
            self._number_of_passes = 0 
            self._f = 0
            self._E = None
            self._BETA  = self.compute_beta()


        @logger.log_decorator
        def compute_beta(self):
            return 1/ (SimulatedTempering.Temperature.k_Boltzmann * self._VALUE)
            # __future__ division -> floting point division


        @logger.log_decorator
        def update_f(self, Tprev) : 
            try:
                self.compute_f(Tprev)
            except SimulatedTempering.Temperature.NoECurrent:
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
                raise SimulatedTempering.Temperature.NoECurrent


        @logger.log_decorator
        def update_E(self, E_new):
            self._number_of_passes += 1 
            try:
                self._E =  self._E  + ( (E_new - self._E ) / self._number_of_passes)
            except TypeError: 
                #self._E = None becase no updated yet
                #First time init : 
                self._E =  ( E_new  / self._number_of_passes)
            
    @logger.log_decorator
    def __init__(self, num_simu, Tmin, Tmax, Tnum, simu_type='md', st_mdp_template_filename = None, **kwargs):
        
        super(SimulatedTempering,self).__init__()
        self._NUM_SIMU = num_simu
        self._T_RANGE=SimulatedTempering.TRange() 
        self._ST_MDP_TEMPLATE_FILENAME= st_mdp_template_filename
        T_range = np.logspace(np.log10(Tmin), np.log10(Tmax), num=Tnum, endpoint=True)
        for T in T_range:
            if simu_type == 'md' : 
                self.create_mdp(T)
            self._T_RANGE.append(SimulatedTempering.Temperature(T))
        for T in self._T_RANGE : 
            print (T._VALUE, T._BETA)
        import pdb; pdb.set_trace()
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
        except SimulatedTempering.TRange.NegativeIndexError : #no previous T because Tcurrent = Tmin 
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
            direction= SimulatedTempering.toss_coin()
            print ('direction = '+str(direction))
            
            T_attempt = self._T_RANGE[self.T_current_idx + direction ]
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
    def write_results (self, idx, E_current_average) : 
        with open('./st.results', 'a') as fout: 
            to_write = "{0}\t{1}\t{2}\t{3}".format(
                idx, 
                self.T_current._VALUE, 
                E_current_average, 
                self.T_current._E)
            for T in self._T_RANGE : 
                to_write += "\t{0}".format(T._f)
            to_write += '\n'
            fout.write (to_write )

    @logger.log_decorator
    def run(self):
        for step_idx in xrange(self._NUM_SIMU) : 
            
            E_current_average = self._SIMULATION.run()
            
            self.T_current.update_E(E_current_average) 

            self.update_f_current()
            self.update_f_next ()

            self.write_results(step_idx, E_current_average)

            T_attempt = self.choose_T_attempt()
        
            if self.attempt_OK(T_attempt) : 
                self._SIMULATION.T_current = T_attempt._VALUE #if MD : change velocity --> Overriding
            else : 
                self._SIMULATION.T_current = self.T_current._VALUE #if MD : change velocity --> Overriding


