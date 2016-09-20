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


```bash
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
                        [--clean-all] [-v]
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
                        md - Molecular Dynamics (currently the only working choice)
                        mc - Monte Carlo (currently not working)
*  --st-mdp-template-filename ST_MDP_TEMPLATE_FILENAME
                        mdp file to use for the ST experiment
*  --st-outname ST_OUTNAME
                        template name for output of the ST experiment
*  --minimisation        Run a minimisation before ST experiment
*  --no-minimisation     Do run a minimisation before ST experiment
*  --minimisation-mdp-filename MINIMISATION_MDP_FILENAME
                        mdp file to use for minimisation
                        _ IF MINIMISATION _
*  --minimisation-outname MINIMISATION_OUTNAME
                        template name for output of minimisation
                        _ IF MINIMISATION _
*  --gene-veloc-outname GENE_VELOC_OUTNAME
                        template name for output of velocities generation
                        _ IF MINIMISATION _
*  --out-path OUT_PATH   Where the outputed results files should be store
                        _/!\ must finish by '/'_
*  --maxwarn MAXWARN     The max number of warnigs allowed when running MD
*  --clean-all           Clean gromacs outputs unecessary to plot the figuresin
                        N'Guyen 2013 /!\ Be carefull
                        _NOT WORKING_
*  -v, --verbose         Turn on detailed info log




Working example : 

```bash 
python ST-on-the-fly.py \
--Tmin 1 --Tmax 5 --Tnum 1 \
--gro-filename ../data/ala10_md000.gro \
--top-filename ../data/ala10.top \
--simu-type md \
--nb-run 3   --st-mdp-template-filename ../data/md1.mdp  --st-outname outst  \
--minimisation \
--minimisation-mdp-filename ../data/mini2.mdp  --minimisation-outname miniout \ 
--maxwarn 1     --out-path ./     -vvv
```


### using Molecular Dynamics with GROMACS

#### without minimisation 


python ST-on-the-fly.py [-h] --Tmin TMIN --Tmax TMAX --Tnum TNUM
                        --gro-filename GRO_FILENAME --top-filename
                        TOP_FILENAME --nb-run NB_RUN --simu-type md
                        --st-mdp-template-filename ST_MDP_TEMPLATE_FILENAME
                        --st-outname ST_OUTNAME
                        --no-minimisation
                        --out-path OUT_PATH [--maxwarn MAXWARN]
                        [--clean-all] [-v]
                        
#### with minimisation 

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
                        [--clean-all] [-v]
                        
                        
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
* xvg file

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

#### Charts : 

TODO : charts



# References

Communication: Simulated tempering with fast on-the-fly weight determination.
Nguyen PH1, Okamoto Y, Derreumaux P. 2013




