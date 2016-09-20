"""
Run GROMACS Molecular Dynamics using a Python wrapper 
"""
#------------------------------------------------------------------------------
# 1. standard library imports
#---
import doctest
import pdb
import shutil 
import shlex
import subprocess

# 2. other imports
#---


# 3. local imports
#---
import logger
#------------------------------------------------------------------------------


class MolecularDynamics(object):
    """MolecularDynamics class 

    Instance Attributs : 
    --------------------
    _mdp_filename : string, private 
        Path to the mdp file 
    _gro_filename : string, private 
        Path to the gro file
    _top_filename : string, private 
        Path to the top file
    _out_path : string, private 
        Path where to store outputs 
    _out_name : string, private 
        Templace for the outputs
    _maxwarn : int, private
        Number of warnings allowed for GROMACS grompp function

    """

    @logger.log_decorator
    def __init__(self, mdp_filename, gro_filename, top_filename, out_path, out_name, maxwarn, **kwargs):
        """ Molecular Dynamics constructor 

        Arguments : 
        -----------
        mdp_filename : string, private 
            Path to the mdp file 
        gro_filename : string, private 
            Path to the gro file
        top_filename : string, private 
            Path to the top file
        out_path : string, private 
            Path where to store outputs 
        out_name : string, private 
            Templace for the outputs
        maxwarn : int, private
            Number of warnings allowed for GROMACS grompp function
        """
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
        """Getter of mdp filename

        Return : 
        --------
        mdp_filename : string
            Path to the mdp file
        """
        return self._mdp_filename
    @property 
    def gro_filename(self):
        """Getter of gro filename

        Return : 
        --------
        gro_filename : string
            Path to the gro file
        """
        return self._gro_filename
    @property 
    def top_filename(self):
        """Getter of top filename

        Return : 
        --------
        top_filename : string
            Path to the top file
        """
        return self._top_filename
    @property 
    def out_path(self):
        """Getter of out path 

        Return : 
        --------
        out_path : string
            Path where to store outputs
        """
        return self._out_path
    @property 
    def out_name(self):
        """Getter of out_name

        Return : 
        --------
        out_name : string
            Template name for the outputs 
        """
        return self._out_name
    @property 
    def maxwarn(self):
        """Getter of maxwarn

        Return : 
        --------
        maxwarn : int
            Number of warning allowed for GROMACS grompp function
        """
        return self._maxwarn


    @logger.log_decorator
    def gmx_grompp(self) : 
        """Run GROMACS grompp command
        """
        cmd = "gmx grompp -f {0} -c {1} -p {2} -o {3}.tpr -po {3}.mdp -maxwarn {4}".format(   \
            self.mdp_filename, \
            self.gro_filename, \
            self.top_filename, \
            self.out_path + self.out_name , \
            self.maxwarn \
        )
        print (cmd)
        pdb.set_trace()
        p = subprocess.Popen(shlex.split(cmd))
        p.wait()
        

    @logger.log_decorator
    def gmx_mdrun(self):
        """Run GROMACS mdrun command
        """
        cmd = "gmx mdrun -v -s {0}.tpr -o {0}.trr  -x {0}.xtc -e {0}.edr -g {0}.log -c {0}.gro".format(self.out_path + self.out_name)
        print (cmd)
        pdb.set_trace()
        p = subprocess.Popen(shlex.split(cmd))
        p.wait()
        

    @logger.log_decorator
    def run(self): 
        """Run a MD using GROMACS grompp then mdrun 
        """
        pdb.set_trace()
        print ('ruuuuuuuuuuuuuun MD')
        self.gmx_grompp()
        self.gmx_mdrun()


################################################################################
if __name__ == "__main__":
    doctest.testmod()