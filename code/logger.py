#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Log module
"""
#------------------------------------------------------------------------------
# 1. standard library imports
#---
import doctest
import functools 
import logging
from logging.handlers import RotatingFileHandler
import time 

# 2. other imports
#--- 

# 3. local imports
#---

#------------------------------------------------------------------------------

# Turn True to have info about function that take 
# the decorator @logger.log_decorator
__DECO_ACTIVATED = False 


__logger = None 


def set_functional_logger(logging_level = logging.DEBUG ):
    """ Set a functional global logger named __logger. 

    Arguments : 
    -----------
    logging_level : logging level - optional
        Threshold for the logger among : 
        logging.CRITICAL -or- 50
        logging.ERROR -or- 40
        logging.WARNING -or- 30
        logging.INFO -or- 20
        logging.DEBUG (default) -or- 10
        logging.NOTSET -or- 0 
    
    Side effect : 
    --------------
    Create a file named 'ST-on-the-fly.log' where logs are appended.
    If this log file reaches 1mo, it gets rolled over.

    Logger details : 
    ----------------
    The logger is formatted as the following :
    time :: levelname :: message

    Example : 
    ----------

    >>> set_functional_logger(logging.INFO)
    >>> __logger.info('there is a logger set to INFO level')
    there is a logger set to INFO level

    """
    global __logger
    __logger = logging.getLogger()
    # set level : 
    __logger.setLevel(logging_level)
    # Add time, level, log message : 
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')

    # File redirection 
    file_handler = RotatingFileHandler( 'ST-on-the-fly.log', 'a', 1000000, 1) 
    #                                                                = 1Mo
    file_handler.setLevel(logging_level)
    file_handler.setFormatter(formatter)
    __logger.addHandler(file_handler)
    
    # Steam (console) redirection : 
    steam_handler = logging.StreamHandler()
    steam_handler.setLevel(logging_level)
    __logger.addHandler(steam_handler)
    



def log_decorator(func):
    """ This decorator prints infos about functions 

    Usage : 
    -------
    Add the decorator @logger.log_decorator before every function on which 
    you want information.
    Turn True the flag __DECO_ACTIVATED

    Side effect : 
    -------------
    Logs module name, function name, args and kwargs at the beginig of the fn
    Logs module name, function name, args, kwargs and timing at the end of the fn

    """
    if not __DECO_ACTIVATED:
        return func
        
    module_name = func.__module__    
    func_name = func.__name__ 
    
    # @functools.wraps reports attributs of the original function on the wrapper
    # (including docstrings)
    @functools.wraps(func) 
    def wrapped(*args, **kwargs):
        msg = "Module={} Function={} Args={} Kwargs={}".format(
            module_name, func_name, args, kwargs
        )       
        __logger.debug("BEGIN " + msg)
        
        t = time.clock()
        data = func(*args, **kwargs)
        
        __logger.debug("END " + msg  + ' Timing : ' + str( time.clock()-t ) )
        return data
    return wrapped



def get_level(v):
    """Get logging level corresponding to verbosity value

    Argument : 
    ----------
    v : int 
        verbosity value

    Return : 
    --------
    level : int 
        logging level 

    Nota Bene : 
    -----------
    The correspondance used is : 
    verbosity | level
    ----------|---------------
        0     |   0 : NOTSET         
        1     |  30 : WARNING
        2     |  20 : INFO
        3     |  10 : DEBUG
    """
    verbosity = [ 0, 30  ,20 ,10]
    try:
        return verbosity[v] 
    except IndexError:
        return 10 


if __name__ == "__main__":
    doctest.testmod()