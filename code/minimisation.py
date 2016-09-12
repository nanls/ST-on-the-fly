#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run minimisation using gromax
"""

from __future__ import print_function

import subprocess
import shlex

import errno
import os


def create_dir_if_not_exists(dir_path): 
    """ Create a dir named <dir_path> if it does not already exists. 
    Parameter
    ---------
    dir_path : str
        Path of the wanted new dir. 
    Side effect
    -----------
    Create a dir named <dir_path> if it does not already exists.
    """

    try:
        os.makedirs(dir_path) # creates all the intermediate directories
    except OSError as e:
        if e.errno == errno.EEXIST :
            print ("The dir already exists.")
        if e.errno == errno.EACCES :
            print ("Permission denied ! ")

if __name__ == "__main__":

    #------
    # run script : 
    cmd = 'bash ./minimisation.sh'
    p = subprocess.Popen(shlex.split(cmd))
    output = p.communicate()
    print (output)

