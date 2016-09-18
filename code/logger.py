



import logging
from logging.handlers import RotatingFileHandler

import time 

import functools 

__DECO_ACTIVATED = True
__logger = None 


def set_functional_logger(logging_level = logging.DEBUG ):
    """ Set a functional global logger named __logger. 

    Arguments
    ---------
    logging_level : logging level - optional
        Threshold for the logger among : 
        logging.CRITICAL
        logging.ERROR
        logging.WARNING
        logging.INFO
        logging.DEBUG (default)
    
    Side effect 
    -----------
    Create a file named 'ST-on-the-fly.log' where logs are appended.

    Example
    -------

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
    if not __DECO_ACTIVATED:
        return func
        
    module_name = func.__module__    
    func_name = func.__name__ 
    
    @functools.wraps(func) # reports attributs of the original function on the wrapper (including docstrings)
    def wrapped(*args, **kwargs):
        msg = "Module={} Function={} Args={} Kwargs={}".format(module_name, func_name, args, kwargs)       
        __logger.debug("BEGIN " + msg)
        
        t = time.clock()
        data = func(*args, **kwargs)
        
        __logger.debug("END " + msg  + ' Timing : ' + str( time.clock()-t ) )
        return data
    return wrapped



def get_level(v):
    verbosity = [ 0, 30  ,20 ,10]
    try:
        return verbosity[v] 
    except IndexError:
        return 10 
