#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Run ST-on-the-fly"""



from __future__ import print_function


class MolecularDynamics(object):
    """docstring for MolecularDynamics"""
    def __init__(self, ):
        super(MolecularDynamics, self).__init__()
        step_duration
        md_duration
        self.T_current


class Simulation(object):
    """docstring for Simulation"""
    def __init__(self):
        super(Simulation,self).__init__()

    def run():
        return compute_E_average()

    def compute_E_average():
        pass

class MonteCarlo(Simulation):
    """docstring for MonteCarlo"""
    def __init__(self, arg):
        super(MonteCarlo, self).__init__()
        self.arg = arg

    def run(): 
        pass


class SimulatedTempering(object):
    """docstring for ST"""

    def __init__(self, ):
        super(SimulatedTempering,self).__init__()
        self._MAX_NUM_STEP
        self._T_RANGE=range(Tmin, Tmax, Tstep)
        self.simulation #pattern strategy
        self.f_current
        self._step_idx
        
        

    def f_attempt_estimate(self, T_attempt):
        # self.simulation.T_current
        # self.simulation.E_average
        pass

    def toss_coin(self):
        pass

    def run_simulation(): #pattern strategy
        self.simulation.run()

    def choose_T_attempt():
        if toss_coin() == 'up' and T_current != Tmax: 
        elif toss_coin() == 'down' and T_current != Tmin : 
            # down 
        else : 
            #do not change T_current

    def compute_metropolis_criterion(T_attempt) : 
        f_attempt = f_attempt_estimate(T_attempt)


        return min (1, ... self.f_current ...  . )

    def attempt_OK():
        mc =  compute_metropolis_criterion() 
        if mc == 1 : 
            return True
        else : 
            return False


    def run(self):
        while self.step_idx < self.MAX_NUM_STEP : 

            run_simulation()

            T_attempt = chose_T_attempt()
        
            if attempt_OK(T_attempt) : 
                self.simulation.T_current = T_next #if MD : change velocity --> Overriding

if __name__ == "__main__":
	print ('go')