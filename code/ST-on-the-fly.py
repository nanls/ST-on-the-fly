#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Run ST-on-the-fly"""



from __future__ import print_function

import logger


import subprocess
import shlex

import random

class MolecularDynamics(object):
    """docstring for MolecularDynamics"""
    @logger.log_decorator
    def __init__(self, mdp_filename, gro_filename, top_filename, out_path, out_name, maxwarn, **kwargs):
        print ('je suis dans le constructor de MolecularDynamics')
        super(MolecularDynamics, self).__init__(**kwargs)
        self._mdp_filename = mdp_filename
        self._gro_filename = gro_filename
        self._top_filename = top_filename
        self._out_path = out_path
        self._out_name = out_name
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
        self.T_current = T_current


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
            print splitted_line
            if splitted_line and splitted_line[0] == 'Potential' : 
                E_average = splitted_line[1]
                return E_average
        logger.__logger.error('E_average could not be found') 
        return 0


@logger.log_decorator
def create_simulation(simu_type, **kwargs): 
    if simu_type == 'md' : 
        return MolecularDynamicsProduction(**kwargs)
    elif simu_type == 'mc':
        return MonteCarlo(**kwargs)
    else : 
        logger.__logger.error('Wrong choice of simulation')
        return None


class SimulatedTempering(object):
    """docstring for ST"""
    @logger.log_decorator
    def __init__(self, num_step, Tmin, Tmax, Tstep, simu_type='md', **kwargs):
        
        super(SimulatedTempering,self).__init__()
        self._NUM_STEP = num_step
        self._T_RANGE=range(Tmin, Tmax, Tstep)
        kwargs['T_current'] = Tmin
        self._SIMULATION=create_simulation(simu_type, **kwargs ) #pattern strategy
        self.f_current=0
        self._step_idx=0
        


    @logger.log_decorator
    def f_attempt_estimate(self, T_attempt):
        # self.simulation.T_current
        # self.simulation.E_average
        pass

    @logger.log_decorator
    @staticmethod
    def toss_coin():
        return random.choice ( [-1, 1] )


    @logger.log_decorator
    def choose_T_attempt(self):
        i_current = self._T_RANGE.index(self._SIMULATION.T_current)
        try:
            T_next = self._T_RANGE[i_current + cls.toss_coin() ]
        except IndexError:
            T_next = self._SIMULATION.T_current
        return T_next

    @logger.log_decorator
    def compute_metropolis_criterion(self, T_attempt) : 
        f_attempt = f_attempt_estimate(T_attempt)

        #min (1, ... self.f_current ...  . )
        return 0
    @logger.log_decorator
    def attempt_OK(self):
        mc =  compute_metropolis_criterion() 
        if mc == 1 : 
            return True
        else : 
            return False

    @logger.log_decorator
    def run(self):
        while self.step_idx < self.MAX_NUM_STEP : 

            self.simulation.run()

            T_attempt = chose_T_attempt()
        
            if attempt_OK(T_attempt) : 
                self.simulation.T_current = T_next #if MD : change velocity --> Overriding

if __name__ == "__main__":
    print ('go')
    logger.set_functional_logger()

    logger.__logger.info('logger OK')

    logger.__logger.info('get args')

    num_step = 100
    dt_pas = 1
    dt_attempt = 10 
    Tmin = 1 
    Tmax = 10 
    Tstep = 1 

    logger.__logger.info('run minimi')
    MD_minimi = MolecularDynamics(        
        mdp_filename = '../data/mini2.mdp', 
        gro_filename = '../data/ala10_md000.pdb', 
        top_filename = '../data/ala10.top', 
        out_path = './', 
        out_name = 'minimisation_before_ST', 
        maxwarn = 13)
    MD_minimi.run()
    logger.__logger.info('minimi OK')

    logger.__logger.info('new ST')
    ST = SimulatedTempering(
        num_step, 
        Tmin, Tmax, Tstep, 
        'md',  
        mdp_filename = '../data/md1.mdp', 
        gro_filename = './minimisation_before_ST.gro', 
        top_filename = '../data/ala10.top', 
        out_path = './', 
        out_name = 'ST-simu', 
        maxwarn = 13
    )
    
    logger.__logger.info('ST.run')


    logger.__logger.info('the end')