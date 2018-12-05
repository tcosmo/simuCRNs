# Welcome!!

SimuCRNs is a simulator for Chemical Reaction Networks (mass action and stochastic). It can build an interactive UI for visualizing CRNs "live".      

# Getting Started

- Run `nix-shell` to prepare the running environment. You might encounters issues on MacOS related to the `send2trash` dependency.    
- Copy `SimuCRNs_demo.ipynb` in the git ignored `bin/` folder.
- Add as a first cell: 
```
import os                         
os.chdir("../")
```
- Have fun :)

# Structure

- `bin/`: git ignored folder for user stuff (notebooks, tests, etc...).    
- `CRNs/`: folder with example CRNs specifications.    
- `default.nix`: nix script to setup the environment.
- `Interaction 101.ipynb`: ipywidget short tutorial, how to build an UI in notebooks.    
- `simuCRNs`: python package for CRN simulation.   
- `SimuCRNs_demo.ipynb`: main demo (c.f. Getting Started).

# TODO

- Adding time in UI    
- Stochastic CRNs
