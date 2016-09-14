#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Run ST-on-the-fly"""



from __future__ import print_function

import logger


class MolecularDynamics(object):
    """docstring for MolecularDynamics"""
    @logger.log_decorator
    def __init__(self, mdp_filename, **kwargs):
        print ('je suis dans le constructor de MolecularDynamics')
        super(MolecularDynamics, self).__init__(**kwargs)
        self._mdp_filename = mdp_filename
    @logger.log_decorator
    def run():
        print ('MD.run()=========')
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
    def run():
        return compute_E_average()
    @logger.log_decorator
    def compute_E_average():
        return 0

class MolecularDynamicsProduction(Simulation,MolecularDynamics):
    """docstring for MolecularDynamicProduction"""
    @logger.log_decorator
    def __init__(self, **kwargs):
        print ('je suis dans le constructor de Prod') 
        super(MolecularDynamicsProduction, self).__init__(**kwargs)

    @logger.log_decorator
    def run():
        super(MolecularDynamicsProduction, self).run() # call MolecularDynamics.run()
        return compute_E_average()
    @logger.log_decorator
    def compute_E_average():
        return 0


@logger.log_decorator
def create_simulation(simu_type, T_current, mdp_filename): 
    if simu_type == 'md' : 
        return MolecularDynamicsProduction(T_current=T_current, mdp_filename=mdp_filename)
    elif simu_type == 'mc':
        return MonteCarlo(T_current)
    else : 
        logger.__logger.error('Wrong choice of simulation')
        return None


class SimulatedTempering(object):
    """docstring for ST"""
    @logger.log_decorator
    def __init__(self, num_step, Tmin, Tmax, Tstep, simu_type='md', mdp_filename=None, **kwargs):
        print (num_step, Tmin, Tmax, Tstep, simu_type, mdp_filename)
        super(SimulatedTempering,self).__init__()
        self._NUM_STEP = num_step
        self._T_RANGE=range(Tmin, Tmax, Tstep)
        self._SIMULATION=create_simulation(simu_type, Tmin, mdp_filename) #pattern strategy
        self.f_current=0
        self._step_idx=0
        


    @logger.log_decorator
    def f_attempt_estimate(self, T_attempt):
        # self.simulation.T_current
        # self.simulation.E_average
        pass
    @logger.log_decorator
    def toss_coin(self):
        pass
    @logger.log_decorator
    def run_simulation(): #pattern strategy
        self.simulation.run()
    @logger.log_decorator
    def choose_T_attempt():
        if toss_coin() == 'up' and T_current != Tmax: 
            print ('up')
        elif toss_coin() == 'down' and T_current != Tmin : 
            print ('down')
            # down 
        else : 
            print ('T does not change')
            #do not change T_current
    @logger.log_decorator
    def compute_metropolis_criterion(T_attempt) : 
        f_attempt = f_attempt_estimate(T_attempt)

        #min (1, ... self.f_current ...  . )
        return 0
    @logger.log_decorator
    def attempt_OK():
        mc =  compute_metropolis_criterion() 
        if mc == 1 : 
            return True
        else : 
            return False

    @logger.log_decorator
    def run(self):
        while self.step_idx < self.MAX_NUM_STEP : 

            run_simulation()

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

    logger.__logger.info('new ST')
    ST = SimulatedTempering(num_step, Tmin, Tmax, Tstep, 'md', 'mdp_file')
    logger.__logger.info('ST.run')


    logger.__logger.info('the end')