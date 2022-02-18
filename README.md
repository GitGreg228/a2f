# a2f - Superconducting properties calculation
Solving Eliashberg equations is effective and accurate way to calculate superconducting critical temperature. The Eliasberg function a2f(w) is required for this method. It can be computed using the parameters of electron phonon coupling. 

The script presented here reads the outputs of the Quantum Espresso with these parameters and computes a2f(w). 

It also calulates:
1. Logarithmic average frequency wlog
2. Mean square frequency w2
3. McMillan and Allen-Dynes superconducting critical temperature Tc
4. Eliashberg superconducting critical temperature Tc (using a2f kernels)
---
## Usage
Run command:
```
python main.py
```
The script will read all files named `*dyn*.elph*` and `output.ph.*` in the directory. From `*dyn*.elph*` files in will take lambas. After that, it will use `output.ph.*` to compute the weights of q-points and invetigate the crystal structure. 

Script parameters are:
1. `-p` - path to the directory with `*dyn*.elph*` and `output.ph.*` files (default: '.')
2. `-s` - exponential smoothing parameter in THz, used to remove acoustic frequancies (default: 3)
3. `-r` - desired resolution of the a2f function (default: cumulative number of positive frequencies in all `*dyn*.elph*` files)
4. `-g` - sigma in gaussian filter used for smoothing (default: 1)
5. `--mu` - Coulomb pseudopotential (default: 0.1)
6. `--tol` - structure tolerance in angstrom (default: 0.2)
