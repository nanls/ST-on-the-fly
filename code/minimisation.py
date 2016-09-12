#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run minimisation using gromax
"""

if __name__ == "__main__":
        
    import subprocess
    import shlex

    cmd = 'bash ./minimisation.sh'
    p = subprocess.Popen(shlex.split(cmd))
    output = p.communicate()
    print (output)

