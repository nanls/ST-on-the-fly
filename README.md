# ST-on-the-fly
Simulated Tempering with fast on-the-fly weight determination. Inspired from Nguyen 2013.

# Disclaimer
This work is still ongoing, and seems to not work properly at the moment.

# Overview
The goal of this project is to provide a tool to run Simulated Tempering.
It implements algorithm described in N'Guyen 2013.

# Contributor :
This work was initiated for a short Python project given in 2nd year of master 
deggree of bioinformatics.
Please contact me at vaginay.athenais@gmail.com for all inquiries.

# Requirements

* [Python2.7](https://www.python.org/download/releases/2.7/) (should works in [Python3](https://www.python.org/download/releases/3.0/) too)

* Python libraries/packages/ecosystems:
  
  * [NumPy](http://www.numpy.org/)
  * [SciPy](http://www.scipy.org/)
  * [Matplotlib](http://matplotlib.org/)
    
  From standard : 
  * [abc](https://docs.python.org/2/library/abc.html)
  * [argparse](https://docs.python.org/2/library/argparse.html)
  * [doctest](https://docs.python.org/2/library/doctest.html)
  * [errno](https://docs.python.org/2/library/errno.html)
  * [functools](https://docs.python.org/2/library/functools.html)
  * [__future__](https://docs.python.org/2/library/__future__.html)
  * [logging](https://docs.python.org/2.7/library/logging.html)
  * [math](https://docs.python.org/2.7/library/math.html)
  * [os](https://docs.python.org/2.7/library/os.html)
  * [pdb](https://docs.python.org/2.7/library/pdb.html)
  * [random](https://docs.python.org/2.7/library/random.html)
  * [re](https://docs.python.org/2/library/re.html)
  * [shlex](https://docs.python.org/2.7/library/shlex.html)
  * [shutil](https://docs.python.org/2.7/library/shutil.html)
  * [subprocess](https://docs.python.org/2.7/library/subprocess.html)
  * [sys](https://docs.python.org/2.7/library/sys.html)
  * [time](https://docs.python.org/2.7/library/time.html)

* To run MD : [GROMAX5.1.2](http://manual.gromacs.org/documentation/5.1.2/)


TODO : check if the code is functional with Python3


# Installation
Just download the archive in <path> and untargz it : 

```bash
cd <path> 
tar zxvf VAGINAY.athenais.tar.gz 
```

# Testing
For testing (doctest) the differnet modules : 

```bash
cd code
python <module_name>
```

TODO : finish unit test

# Usage

## Run Simulated Tempering 


### using Molecular Dynamics with GROMACS

### using Monte Carlo

TODO

## Run Molecular Dynamics

TODO 

## Run Monte Carlo experiments

TODO

## Output visualisation

TODO




# File formats

## Input : 

### Classical GROMACS files

* gro file
* top file
* mdp file

You can find a description of the format [here](http://manual.gromacs.org/online/files.html)
## Output  :

### Output of `code/md.py`

#### Classical GROMACS files

* gro file
* top file
* mdp file
* edr file
* log file
* tpr file
* xtc file

You can find a description of the format [here](http://manual.gromacs.org/online/files.html)

### Output of `code/st.py`


#### Classical GROMACS files

* gro file
* top file
* mdp file
* edr file
* log file
* tpr file
* xtc file

You can find a description of the format [here](http://manual.gromacs.org/online/files.html)

#### Results file : 

* .results

tab separated file with following fields : 

* idx : int < <num-run> - attempt index 
* t_current : float - time at the begining of the run 
* T_current : float - temperature of the run 
* E_run : float - energy of the run 
* E_T : float - updated energy for T_current 
* For each temperature in provided T range : f(T) 

Charts



# References




