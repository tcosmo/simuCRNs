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

# TODO

- Adding time in UI    
- Stochastic CRNs
