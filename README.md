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
Please contact Athenais VAGINAY at vaginay.athenais@gmail.com for all inquiries.

# Requirements

* [Python2.7](https://www.python.org/download/releases/2.7/) 
    (should works in [Python3](https://www.python.org/download/releases/3.0/) too)

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

Each module can be tested (doctest) running the module as a script : 

```bash
cd code
python <module_name>
```

This causes the examples in the docsttrings to get executed and verified. 

TODO : finish unit test

# Code info: 

The code is not fully PEP8 compliant.

## Documentation  :

You can use pydoc outputed documentation. 
html files are store in `./code/doc`

# Usage


/!\
Not not cheat on arguments because combinations are not verified in the code.


TODO : verify arguments combination 


The script should be able to recognize abbreviation of long option.

## Run Simulated Tempering (ST)



```bash
cd code
python ST-on-the-fly.py [-h] --Tmin TMIN --Tmax TMAX --Tnum TNUM
                        [--gro-filename GRO_FILENAME] [--top-filename
                        TOP_FILENAME] --nb-run NB_RUN --simu-type {md,mc}
                        [--st-mdp-template-filename ST_MDP_TEMPLATE_FILENAME]
                        --st-outname ST_OUTNAME
                        [--minimisation | --no-minimisation]
                        [--minimisation-mdp-filename MINIMISATION_MDP_FILENAME]
                        [--minimisation-outname MINIMISATION_OUTNAME]
                        [--gene-veloc-outname GENE_VELOC_OUTNAME]
                        --out-path OUT_PATH [--maxwarn MAXWARN]
                        [--clean-all] [-v[vv]]
```


Arguments Explanation :

*  -h, --help            show this help message and exit
*  --Tmin TMIN           At which temperature (in Kelvin) begins the experiment
*  --Tmax TMAX           The maximal temperature (in Kelvin) during the
                        experiment
*  --Tnum TNUM           Number of temperature in the given range [Tmin, Tmax]
*  --gro-filename GRO_FILENAME : 
                        _gro_ file to use for the ST experiment
*  --top-filename TOP_FILENAME : 
                        topology file to use for the ST experiment
*  --nb-run NB_RUN       The number of molecular dynamics during the experiment
                        :t_one-md * num-run = t_ST
*  --simu-type {md,mc}   type of simulation : 
                        md - Molecular Dynamics 
                            (currently the only working choice)
                        mc - Monte Carlo 
                            (currently not working)
*  --st-mdp-template-filename ST_MDP_TEMPLATE_FILENAME
                        mdp file to use for the ST experiment
*  --st-outname ST_OUTNAME
                        template name for output of the ST experiment
                        Should not match either minimisation-outname 
                        nor gene-veloc-outname
*  --minimisation        Run a minimisation before ST experiment
*  --no-minimisation     Do run a minimisation before ST experiment
*  --minimisation-mdp-filename MINIMISATION_MDP_FILENAME
                        mdp file to use for minimisation
                        _ IF MINIMISATION _
*  --minimisation-outname MINIMISATION_OUTNAME
                        template name for output of minimisation
                        Should not match either st-outname nor gene-veloc-outname
                        _ IF MINIMISATION _
*  --gene-veloc-outname GENE_VELOC_OUTNAME
                        template name for output of velocities generation
                        Should not match either st-outname nor minimisation-outname
                        _ IF MINIMISATION _
*  --out-path OUT_PATH   Where the outputed results files should be store
                        _/!\ must finish by '/'_
*  --maxwarn MAXWARN     The max number of warnings allowed when running MD
*  --clean-all           Clean gromacs outputs unecessary to plot the figures in
                        N'Guyen 2013 /!\ Be carefull : _NOT WORKING_
*  -v, --verbose         Turn on detailed info log
                        The level of verbosity is function of the v count






### using Molecular Dynamics with GROMACS

#### without minimisation 

```bash
cd code
python ST-on-the-fly.py [-h] --Tmin TMIN --Tmax TMAX --Tnum TNUM
                        --gro-filename GRO_FILENAME --top-filename
                        TOP_FILENAME --nb-run NB_RUN --simu-type md
                        --st-mdp-template-filename ST_MDP_TEMPLATE_FILENAME
                        --st-outname ST_OUTNAME
                        [--no-minimisation]
                        --out-path OUT_PATH [--maxwarn MAXWARN]
                        [--clean-all] [-v[vv]]
``` 


Working example : 

```bash 
cd code
python ST-on-the-fly.py \
--Tmin 1 --Tmax 5 --Tnum 1 \
--gro-filename ../data/ala10_md000.gro \
--top-filename ../data/ala10.top \
--nb-run 3  --simu-type md \
--st-mdp-template-filename ../data/md1.mdp  --st-outname outst  \
--out-path ../result/test-without-minimi/   --maxwarn 1       -vvv
```

                    
#### with minimisation 


```bash
cd code
python ST-on-the-fly.py [-h] --Tmin TMIN --Tmax TMAX --Tnum TNUM
                        --gro-filename GRO_FILENAME --top-filename
                        TOP_FILENAME --nb-run NB_RUN --simu-type md
                        --st-mdp-template-filename ST_MDP_TEMPLATE_FILENAME
                        --st-outname ST_OUTNAME
                        --minimisation
                        --minimisation-mdp-filename MINIMISATION_MDP_FILENAME
                        --minimisation-outname MINIMISATION_OUTNAME
                        --gene-veloc-outname GENE_VELOC_OUTNAME
                        --out-path OUT_PATH [--maxwarn MAXWARN]
                        [--clean-all] [-v[vv]]
```    

Working example : 

```bash 
cd code
python ST-on-the-fly.py \
--Tmin 1 --Tmax 5 --Tnum 1 \
--gro-filename ../data/ala10_md000.gro \
--top-filename ../data/ala10.top \
--nb-run 3  --simu-type md \
--st-mdp-template-filename ../data/md1.mdp  --st-outname outst  \
--minimisation \
--minimisation-mdp-filename ../data/mini2.mdp  --minimisation-outname miniout \ 
--gene-veloc-outname gene-vel \
--out-path ../result/test-with-minimi/  --maxwarn 1  -vvv 

```     
                        
### using Monte Carlo

Not working 

TODO


## Output visualisation

Use : 

```bash 
python plot.py ... 
```

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

* cpt file
* edr file
* gro file
* log file
* mdp file
* top file
* tpr file
* xtc file

You can find a description of the format [here](http://manual.gromacs.org/online/files.html)

### Output of `code/st.py`


#### Classical GROMACS files


* cpt file 
* edr file
* gro file
* log file
* mdp file
* top file
* tpr file
* xtc file
* xvg file

You can find a description of the format [here](http://manual.gromacs.org/online/files.html)

#### Results file : 

* <out-path>/<st-out-name>.results

tab separated file with following fields : 

* idx : int < <num-run> - attempt index 
* t_current : float - time at the begining of the run 
* T_current : float - temperature of the run 
* E_run : float - energy of the run 
* E_T : float - updated energy for T_current 
* For each temperature in provided T range : f(T) 

#### Charts : 

TODO : charts

# Known Bugs : 

* Probabilities of attempt given by the Metropolis Criteriion is always 0 or 1
could be because energies are to big

* The code is vulnerable to TOCTTOU (Time Of Check to Time Of Use) Error



# References

Communication: Simulated tempering with fast on-the-fly weight determination.
Nguyen PH1, Okamoto Y, Derreumaux P. 2013




