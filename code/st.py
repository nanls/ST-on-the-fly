#!/usr/bin/env python
# -*- coding: utf-8 -*-


#------------------------------------------------------------------------------

# 1. standard library imports
#---

# from __future__ imports must occur at the beginning of the file.
# if not : SyntaxError 
from __future__ import print_function
from __future__ import division 

import abc
import doctest
import math
import os
import pdb
import random
import re
import shlex
import shutil
import subprocess
import sys

# 2. other imports
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


# 3. local imports
#---
import logger
logger.set_functional_logger()
global log
log = logger.__logger

from md import MolecularDynamics

#------------------------------------------------------------------------------

@logger.log_decorator
def create_simulation(simu_type, **kwargs): 
    """ Create a simulation of the chosen type

    Argumuents : 
    ------------
    simu_type : string
        Type of the simulation. 
        Could be : 
        'md' - Molecular Dynamics
        'mc' - Monte Carlo 
        Currently the only one that works is 'md'
    **kwargs :
        Every named argumnents needed for the constructor of the chosen type.

    Return : 
    --------
    None if simu_type does not correspond to a rigth choice.
    -or-
    simulation : Simulation derived object
        Could be :  
        MolecularDynamicsProduction object 
        MonteCarlo object 
        Currently : can only return MolecularDynamicsProduction
    """
    if simu_type == 'md' : 
        return MolecularDynamicsProduction(**kwargs)
    elif simu_type == 'mc':
        return MonteCarlo(**kwargs)
    else : 
        log.error('Wrong choice of simulation')
        return None

class Simulation(object):
    """Simulation class
    
    Corresponds to the first step of ST experiment 

    Instance Attributs : 
    ---------------------
    T_current : int 
        Current temperature

    Must be overrided by child classes: 
    -----------------------------------
    run method
    get_simu_step method
    """
    @logger.log_decorator
    def __init__(self, T_current, **kwargs):
        """Simulation constructor 

        Argument : 
        ----------
        T_current : float 
            current temperature of the simulation 
            (it should be Tmin of the ST experiment)
        """
        print ('je suis dans le constructor de Simulation')

        super(Simulation,self).__init__(**kwargs)
        print ('je suis revenu de super')
        self._T_current = T_current

    @property
    def T_current(self):
        """ Get T_current attribut

        Return : 
        T_current attribut : float 
            current temperature of the simulation 
        """
        return self._T_current

    @T_current.setter
    def T(self, new_T) : 
        """ Set T

        Argument : 
        ----------
        new_T : float > 0 
            New value for T_current 

        Raise : 
        -------
        ValueError if new_T < 0
        """
        if new_T > 0 : 
            self._T_current = new_T
        else : 
            raise ValueError

    @abc.abstractmethod
    def run(self) : 
        """ Run the simulation

        This method has to be override by any derived class
        and should return Energy of the simulation.

        """
        pass
    @abc.abstractmethod
    def get_simu_step(self) : 
        """Get the simulation duration 

        This method has to be overided by any derived class 
        and should return one simulation duration
        """
        pass
    

class MonteCarlo(Simulation):
    """MonteCarlo class 

    A possible kind of first step for ST experiment.
    Currently : not implemented

    Instance Attribut : 
    -------------------
    Those of simulation 

    """
    
    @logger.log_decorator
    def run(self):
        """ Run a Monte Carlo simulation 

        override Simulation

        Return : 
        ---------
        E : float 
            energy of the last Monte Carlo simulation 
        """
        return compute_E_average()
    @logger.log_decorator
    def compute_E_average(self):
        """ Compute E on the last Monte Carlo simulation 

        Return : 
        --------
        E : float 
            energy of the last Monte Carlo simulation 
        """
        return 0 # ??? XXX pass

class MolecularDynamicsProduction(Simulation,MolecularDynamics):
    """MolecularDynamicProduction class

    A possible kind of first step for ST experiment.
    Currently: working using GROMACS 

    Instance Attribut : 
    -------------------
    Those of Simulation 
    Those of MolecularDynamics

    MDP_TEMPLATE : string , constant
        Path to the mdp template file 

    """
    @logger.log_decorator
    def __init__(self, T_range, **kwargs):
        """constructot of MolecularDynamicsProduction instance

        Arguments : 
        ----------
        T_range : list of float 
            list of temperature value of the ST experiment 
        kwargs : dict
            every arguments needed for parent constructors
        """
        print ('je suis dans le constructor de Prod') 
        super(MolecularDynamicsProduction, self).__init__(**kwargs)
        self._MDP_TEMPLATE = kwargs['mdp_filename'] 
        for T in T_range : 
            self.create_mdp(T)
        self.mdp_filename = '{0}_{1}.mdp'.format(self.MDP_TEMPLATE , self.T_current)

    @property
    def MDP_TEMPLATE(self):
        return self._MDP_TEMPLATE

    

    @logger.log_decorator
    def create_mdp(self, T) : 
        """Create an mdp file for the given Temperature 

        Arguments : 
        ----------
        T : float 
            Temperature value 

        Side effect : 
        -------------
        Create a file named <args.st_mdp_template_filename>_<T>.mdp
        where the field gen_temp if set to <T>
        """
        
        mdp_filename = '{0}_{1}.mdp'.format(self.MDP_TEMPLATE, T)

        with open(self.MDP_TEMPLATE, 'r') as infile, \
            open(mdp_filename, 'w') as outfile    :
            for line in infile : 
                if line.startswith('gen_temp') : 
                    line = "gen_temp = {0}\n".format(T)

                outfile.write(line)



    @logger.log_decorator
    def get_simu_step(self):
        """Get the duration of an MDP looking in mdp template file

        overide Simulation 

        Return : 
        ---------
        duration : float 
            duration of an MDP
        """
        with open(self.MDP_TEMPLATE, 'r') as fin: 
            for line in fin : 
                if line.startswith("dt") : 
                    dt= float(re.split(r'[=;]', line)[1])
                elif line.startswith("nsteps") : 
                    nsteps = float(re.split(r'[=;]', line)[1])
        return dt * nsteps

    @logger.log_decorator
    def run(self, tcurrent):
        """ Run an MDP 

        Agument : 
        ---------
        tcurrent : float 
            ST time at the begining of the MDP 

        Return : 
        --------
        E : float 
            energy of the run 
        """
        save_out_name = self.out_name
        print (save_out_name)
        self.out_name += str(tcurrent)
        print (self.out_name)
        print ('j apelle super MDrun ')
        MolecularDynamics.run(self) # call MolecularDynamics.run()
        E = self.compute_E_average()

        os.rename("{0}/{1}.gro".format(self.out_path, self.out_name), self.gro_filename)
        self.cat_edr(tcurrent)
        self.cat_xtc(tcurrent)
        self.out_name = save_out_name
        return E


    def cat_gmx_files(self, fn, ext, t_current):
        """Concatenate Gromacs file 

        Arguments : 
        -----------
        fn : string 
            The Gromacs function that has to be used 
        ext : string 
            extention of the file to concatenante 
        t_current : float 
            ST time at the begining of the MDP 

        Side effect :
        -------------
        Concatenate desired files within <out_path>/cat.<out_name>.<ext>
        """
        
        if t_current == 0 : 
            cmd = 'gmx {0} -f {1}/{2}.{3} -o {1}/cat.{3}'.format(fn, self.out_path, self.out_name, ext)
            p = subprocess.Popen(shlex.split(cmd))
            p.wait()
        else : 
            cmd1 = "echo '0\n{0}\n'".format(t_current)
            p1 =  subprocess.Popen(cmd1, shell = True, stdout=subprocess.PIPE)

            cmd2 = ("gmx {0} -f {1}/cat.{3} {1}/{2}.{3} -o {1}/cat.{3}  -settime ".format(fn, self.out_path, self.out_name, ext))

            p2 = subprocess.Popen(shlex.split(cmd2), stdin=p1.stdout)
            p2.wait()
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        
       
    def cat_xtc(self, t_current):
        """Concatenate xtc files using trjcat then rm the file of the run 

        Argument : 
        ---------
        t_current : float 
            ST time at the begining of the MDP 

        """
        self.cat_gmx_files('trjcat', 'xtc', t_current)
        os.remove("{0}{1}.xtc".format(self.out_path, self.out_name))

    def cat_edr(self, t_current):
        """Concatenate edr files using eneconv then rm the file of the run 

        Argument : 
        ----------
        t_current : float 
            ST time at the begining of the MDP 
        """
        self.cat_gmx_files('eneconv', 'edr', t_current)
        os.remove("{0}{1}.edr".format(self.out_path, self.out_name))


    @logger.log_decorator
    def gmx_energy(self, arg = 'Potential') : 
        """ Run gmx energy on the last run 

        Argumnent : 
        -----------
        arg : string (default : Potential)
            argument given to gmx energy
        """

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
        """Compute E average of the last MDP run using gmx energy 

        Return : 
        --------
        E : float 
            energy of the last MDP run 
        """
        output = self.gmx_energy('Potential')
        for line in output : 
            splitted_line = line.split()
            print (splitted_line)
            if splitted_line and splitted_line[0] == 'Potential' : 
                E_average = splitted_line[1]
                return -int(E_average) 
                # in NGuyen 2013, Energies are defined as integral -> positive
                # in gmx, E are negative -> take the oposite to have positive value
        log.error('E_average could not be found') 
        

    @Simulation.T_current.setter
    def T_current(self, T_new):
        """Set T_current attribut 
        
        Overdide Simulation 

        Argument : 
        ----------
        T_new : float 
            new value of T_current 

        Side effect : 
        -------------
        change velocities in the gro file
        change the mdp file to use 
        """
        
        print ("setter de Tcurrent ============")
        self.update_velocities(T_new)

        self._T_current = T_new 
        # not self.T_current = T_new
        # otherwise it calls the setter that calls the setter, that c... 
        # do not use setter in setter ! 
        
        self.mdp_filename = self.MDP_TEMPLATE % self.T_current

        
        

    def update_velocities(self, T_new):
        """Update velocities in the gro file 

        Argument : 
        ----------
        T_new : float 
            new temperature value 

        Side effect : 
        -------------
        Change the velocities in the gro file 
        Use a tmp file 
        Then replace the old file by the tmp one. 
        """

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
    
    Currently : 
    works only using a Molecular Dynamics 
    Monte Carlo usage is an ongoing work

    Nested classes : 
    ----------------
    TRange
    Temperature 

    Instance Attributes : 
    ---------------------

        NUM_SIMU : int 
            The number of simulation (therefore number of attempt) to do
        T_RANGE : TRange (list-like)
            List of Temperature used in the ST experiment 
        _ST_MDP_TEMPLATE_FILENAME : string 
            Path to the mdp file used as template for the ST experiment 
        simu_step : float 
            Duration of a Molecular Dynamics 
        _SIMULATION : Simulation derived object
            An instance of the simulation of the desired type
        self._RES_FILENAME : string 
            Path to the file where results will be saved

    Static Class Method : 
    ---------------------
    toss_coin
       
    """

    class TRange(list):
        """ TRange class 

        A list-like without negative indiciation

        Specific Exception : 
        --------------------
        NegativeIndexError 

        Attributes : 
        -----------
        Those of list
        """

        class NegativeIndexError(IndexError):
            """NegativeIndexError Exeption

            Atributes : 
            ------------
            Those of IndexError
            """
            def __init__(self,*args,**kwargs):
                super(Exception, self).__init__(*args,**kwargs)


        def __getitem__(self, idx):
            """Get item corresponding to the given index

            Override list getitem to disable negative indices

            Argument : 
            ----------
            idx : int positive 
                index of the wanted element 
                0 <= idx <= len (list)

            Return : 
            --------
            wanted element 

            Raise : 
            -------
            NegativeIndexError : 
                if idx < 0 
            IndexError : 
                if idx > len (list)
            TypeError : 
                if idx not int 
            """
            if isinstance(idx, int):
                if idx < 0:
                    
                    raise SimulatedTempering.TRange.NegativeIndexError ("neg index not allowed")
                    
                else :  
                    return super(SimulatedTempering.TRange, self).__getitem__(idx)
            else : 
                raise TypeError ("TRange indices must be integers, not "+ str(type (idx) ) )  




    class Temperature(object):
        """Temperature class

        Specific Exeption : 
        -------------------
        NoECurrent 

        Class attribute: 
        ----------------
        k_Boltzmann : float 
            Boltzmann constant in kJ / K (same unit than Gromacs)

        Instance attributes : 
        ---------------------
        VALUE : float, constant
            The temperature itself in Kelvin (same unit than Gromacs)
        number_of_passes : int
            The number of time ST used this temperature
        f : float
            The weight associated to this temperature 
        E : float / None
            The average potential energy for this temperature 
        BETA : float, constant 
            the beta constant for this temperature
            = 1 / (k_Boltzmann * VALUE)
        """

        class NoECurrent(Exception):
            """NoECurrentExeption

            Raised when trying to use E that is None 
            """
            def __init__(self):
                super(SimulatedTempering.Temperature.NoECurrent, self).__init__()


        k_Boltzmann = constants.value(u'Boltzmann constant') / 1000 
        # using gromacs, E are in kilo Joule / K 
        # scipy gives Boltmann constant in Joule per K -> / 1000 

        def __init__(self, value):
            """ Temperature constructor

            Argument : 
            ----------
            value : float 
                The value of the temperature
            """
            super(SimulatedTempering.Temperature, self).__init__()
            self._VALUE = value
            self._number_of_passes = 0 
            self._f = 0
            self._E = None
            self._BETA  = self.compute_beta()

        @property
        def VALUE(self):
            """get VALUE

            Return : 
            --------
            VALUE : float 
                Value stored for this Temperature
            """
            return self._VALUE
        @property
        def number_of_passes(self):
            """get number_of_passes

            Return : 
            --------
            number_of_passes : int 
                Counter of how many simulation was run with this temperature value
            """
            return self._number_of_passes
        @number_of_passes.setter
        def number_of_passes(self, new) : 
            """Set number_of_passes

            Argument : 
            ----------
            new : int > 0 
                new value for number_of_passes 

            Raise : 
            -------
            ValueError if new < 0
            """
            if new >= 0 : 
                self._number_of_passes = new
            else : 
                raise ValueError

        @property
        def f(self):
            """ get f 

            Return : 
            ---------
            f : float 
                Weight associated with this temperature
            """
            return self._f
        @f.setter
        def f(self, new_f):
            """ set f 

            Argument : 
            ---------
            new_f : float 
                new f value 
            """
            self._f = new_f
        
        @property
        def E(self):
            """Get E

            Return : 
            --------
            E : float 
                The average potential energy for this temperature 

            Raise :
            -------
            NoECurrent : 
                if E is None

            """
            if not self._E : 
                raise SimulatedTempering.Temperature.NoECurrent
            else : 
                return self._E
        @E.setter 
        def E (self, new): 
            """Set E 

            Argument : 
            ----------
            new : float 
                new E value 
            """
            self._E = new
        @property
        def BETA(self):
            """ Get BETA 

            Return 
            -------
            BETA : float 
                beta constant associated to this temperature
            """
            return self._BETA
        
        
        @logger.log_decorator
        def compute_beta(self):
            """Compute the beta constant of the temperature

            Return : 
            --------
            beta : float 
                beta constant of the Temperature 
            """
            return 1/ (SimulatedTempering.Temperature.k_Boltzmann * self.VALUE)
            # __future__ division -> floting point division


        @logger.log_decorator
        def update_f(self, Tprev) : 
            """Update the weight of the Temperature 

            The update is done with either a computation 
            or with an estimation if the computation is impossible.
            """
            try:
                self.compute_f(Tprev)
            except SimulatedTempering.Temperature.NoECurrent:
                self.estimate_f(Tprev)


        @logger.log_decorator
        def estimate_f(self, Tprev):
            """Do an estimation of the weight of the temperature

            Use the formula in Nguyen 2013 : 
            (beta_next - beta_previous) E_previous / 2
            """

            self.f = (self.BETA - Tprev.BETA ) * Tprev.E / 2



        
        @logger.log_decorator
        def compute_f(self, Tprev):
            """Do a computation of the weight of the temperature 

            Use the formula in Nguyen2013 
            f_prev + (beta_curr - beta_prev) (E_curr + E_prev) / 2 

            Raise : 
            -------
            NoECurrent : 
                If E_Tcur is not available

            Nota Bene : 
            -----------
            Exits if Tprev do not have Energy
            """
            try:
                self.f =  Tprev.f + (self.BETA - Tprev.BETA ) * ( self.E + Tprev.E )  / 2
            except SimulatedTempering.Temperature.NoECurrent:
                if not self.E : 
                    print ('no E for T = {}'.format(self.VALUE))
                    raise 
                elif not Tprev.E : 
                    print ('problem when trying to access E for T = {}'.format(self.VALUE))
                    exit()

        def update_number_of_passes(self):
            """ Increment by one the number of passes
            """
            self.number_of_passes += 1 

        @logger.log_decorator
        def update_E(self, E_new):
            """Update E including E_new in E

            Parameter : 
            ------------
            E_new : float 
                The average potential energy for the last MD runed with this temperature 
                to be included in E.

            """
            try:
                self.E =  self.E  + ( (E_new - self.E ) / self.number_of_passes)
            except SimulatedTempering.Temperature.NoECurrent: 
                #self._E = None becase no updated yet
                #First time init : 
                self.E =  ( E_new  / self.number_of_passes)
            
    @logger.log_decorator
    def __init__(self, num_simu, Tmin, Tmax, Tnum, res_filename, simu_type='md', **kwargs):
        """ Constructor of an instance of SimulatedTempering 

        Aguments : 
        ----------
        num_simu : int 
            The number of simulation (therefore number of attempt) to do
        Tmin : float 
            Value of the temperature the ST experiment begins with 
        Tmax : float 
            Maximal value of temperature the ST experiment can reach
        Tnum : int 
            Number of temperature in the given range [Tmin ; Tmax] 
        res_filename : string
            Path to the file results have to be saved 
        simu_type : string 
            Type of simulaiton to use. 
            Can be : 
            'md' - Molecular Dynamics (by default)
            'mc' - Monte Carlo 
            Currently : 
            'md' is the only working choice.
            
        """
        super(SimulatedTempering,self).__init__()
        self._NUM_SIMU = num_simu
        self._T_RANGE=SimulatedTempering.TRange() 

        T_range = np.logspace(np.log10(Tmin), np.log10(Tmax), num=Tnum, endpoint=True)
        for T in T_range:
            self._T_RANGE.append(SimulatedTempering.Temperature(T))
        for T in self._T_RANGE : 
            print (T.VALUE, T.BETA)
        #range (a, b) = [a, b[
        #range (a, b+1) = [a, b+1[ = [a, b]
        kwargs['T_current'] = self._T_RANGE[0].VALUE
        self._SIMULATION=create_simulation(simu_type, T_range = [T.VALUE for T in self._T_RANGE], **kwargs ) #pattern strategy
        self.simu_step = self._SIMULATION.get_simu_step()

        #T in T_range must be initialized before init res file 
        self._RES_FILENAME =res_filename
        self.init_res_file()

    def init_res_file(self):
        """Init the res file 

        Side effect : 
        -------------
        Create a file named <self._RES_FILENAME>
        /!\ if a file with the same name already exists, it override it
        and write the header 
        """
        with open(self._RES_FILENAME, 'w') as fout : 
            to_write = "idx\tt_current\tT_current\tE_MD\tE_T"
            for T in self._T_RANGE : 
                to_write += "\tf({0})".format(T.VALUE)
                print (to_write)
            to_write += '\n'
            fout.write (to_write )

    @property
    def T_current(self):
        """Getter of T_current 

        Return : 
        --------
        T : Temperature object 
            Temperature object corresponding to the current T
        """
        return self._T_RANGE[self.T_current_idx]
    
    @property    
    def T_current_idx(self):
        """Getter of T_current_idx

        Return : 
        --------
        Tidx : int 
            idx in T_RANGE of the Temperature object corresponding to the current T
        """
        return self.get_T_idx(self._SIMULATION.T_current)

    def get_T_idx(self, T_wanted) : 
        """Get idx in T_RANGE of the wanted temperature 
        
        Argument : 
        ----------
        T_wanted : float 
            Value of the wanted temperature 

        Return : 
        --------
        i : int
            index in T_RANGE of the Temperature object which the value is T

        """
        return [T.VALUE for T in self._T_RANGE].index(T_wanted)
        
    @logger.log_decorator
    def update_f_current(self):
        """ Update weight of the current Temperature object

        Nota Bene : 
        -----------
        f(Tmin) is always equal to 0 
        """
        try : 
            T_previous = self._T_RANGE[self.T_current_idx-1]
            self.T_current.update_f(T_previous)
        except SimulatedTempering.TRange.NegativeIndexError : #no previous T because Tcurrent = Tmin 
            self.T_current.f = 0 #f_Tmin is always equal to 0.

        # Remember : 
        # If an exeption occurs during execution of the try clause,
        # the rest of the clause is skipped.    

    @logger.log_decorator
    def update_f_next(self):
        """Update weight of the upper numbering Temperature object

        Nota Bene : 
        -----------
        if T_current = Tmax, there is not update to do 
        """
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
        """Toss a coin 

        Return : 
        -------
        result : int 
            random choice of -1 or 1
        """
        return random.choice ( [-1, 1] )


    @logger.log_decorator
    def choose_T_attempt(self):
        """Choose T attempt 

        Return : 
        --------
        T_attempt : Temperature object  / None 
            The temperature ST will attempt
            None if no attempt 
        """
        try:
            direction= SimulatedTempering.toss_coin()
            print ('direction = '+str(direction))
            
            T_attempt = self._T_RANGE[self.T_current_idx + direction ]
        except IndexError:
            T_attempt = None
        return T_attempt

    @logger.log_decorator
    def compute_metropolis_criterion(self, T_attempt) : 
        """Compute metropolis criterion 

        Argument : 
        ----------
        T_attempt : Temperature object 
            The temperature of which ST is trying 

        Return : 
        --------
        mc : float 
            Metropolis criterion computed using the formula in Nguyen 2013
            it corresponds to the probability of the choice to move to <T_attempt>
        """

        try:
            insider = - (                                                        \
                ( T_attempt.BETA- self.T_current.BETA )  * self.T_current.E   \
                - ( T_attempt.f - self.T_current.f )                           \
            )           
            res = math.exp ( insider  ) 
        except OverflowError: 
            if insider < 0 : 
                res = 0
            elif insider > 0 : 
                res = float('inf')
            else : 
                print ('case not handled !! ')
        
        mc = min (1, res)
        return mc
        # mc = min (1 , exp (-   [  (beta_attempt - beta_current) * E_current_average  -   (f_attempt_estimate - f_current_compute)  ]   ) )

    @logger.log_decorator
    def attempt_OK(self, T_attempt):
        """Determine wheather the simulation should move to T_attempt or not 

        Argument : 
        ----------
        T_attempt : Temperature object 
            The temperature ST is trying 
        Return : 
        --------
        answer : boolean 
            wheather the simulation will move to T_attempt or not 
        """
        if not T_attempt : return False
        mc =  self.compute_metropolis_criterion(T_attempt) 

        return np.random.choice([True, False], 1, p= [ mc , 1 - mc  ])[0]


    @logger.log_decorator
    def append_results (self, step_idx, E_MD) : 
        """Update the result files with new info about the run

        Arguments: 
        -----------
        step_idx : int 
            run idx
        E_MD : float 
            Average Epot of the run 
        """
        with open(self._RES_FILENAME, 'a') as fout: 
            to_write = "{0}\t{1}\t{2}\t{3}\t{4}".format(
                step_idx, 
                step_idx * self.simu_step, #t_current
                self.T_current._VALUE, 
                E_MD, 
                self.T_current.E)
            for T in self._T_RANGE : 
                to_write += "\t{0}".format(T.f)
            to_write += '\n'
            fout.write (to_write )


    @logger.log_decorator
    def run(self):
        """Run the ST experiment 

        Algorithm : 
        ----------
        The algorithm is fully described in Nguyen 2013. 

        For short : 

        for each run : 
            run a simulation at Tcurrent
            accumulate potential energy for Tcurrent
            estimate f(Tcurrent)
            estimate f(Tnext) 
            choose Tattempt 
            Attempt switching to temperature Tattempt 
            if attempt successful : 
                Tcurrent = Tattempt

        nota bene : the ST starts at the lowest temperature

        """

        for step_idx in xrange(self._NUM_SIMU) : 
            t_current = step_idx * self.simu_step

            E_current_average = self._SIMULATION.run(t_current)
            self.T_current.update_number_of_passes()
            self.T_current.update_E(E_current_average) 

            self.update_f_current()
            self.update_f_next ()

            self.append_results(step_idx, E_current_average)

            T_attempt = self.choose_T_attempt()
        
            if self.attempt_OK(T_attempt) : 
                self._SIMULATION.T_current = T_attempt.VALUE #if MD : change velocity --> Overriding


################################################################################
if __name__ == "__main__":
    doctest.testmod()