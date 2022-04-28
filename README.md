# a2f - Superconducting properties calculation
Solving Eliashberg equations is effective and accurate way to calculate superconducting critical temperature. The Eliasberg function a2f(w) is required for this method. It can be computed using the parameters of electron phonon coupling. 

The script presented here reads the outputs of the Quantum Espresso with these parameters and computes a2f(w). 

It also calculates:
1. Logarithmic average frequency wlog
2. Mean square frequency w2
3. McMillan and Allen-Dynes superconducting critical temperature Tc
4. Eliashberg superconducting critical temperature Tc (using a2f kernels)
5. Parameters of superconducting state:
6.1 
---
## Usage
In `src/` dir, run command:
```
python main.py -p ../desired/path
```
The script will read all files named `*.dyn*.elph*` and `output.ph.*` in the `../desired/path` directory. First, lambdas and gammas are taken from `*.dyn*.elph*`. After that, `output.ph.*` will be used to compute the weights of q-points and investigate the crystal structure. 

For more detailed explanation, see [notebooks/a2f_tutorial.ipynb](Tutorial)

Script input parameters are:
1. `-p` - path to the directory with `*.dyn*.elph*` and `output.ph.*` files (default: `.`)
2. `-s` - exponential smoothing parameter in THz, used to remove acoustic frequencies (default: 3)
3. `-r` - desired resolution of the a2f function (default: cumulative number of positive frequencies in all `*dyn*.elph*` files)
4. `-g` - sigma in gaussian filter used for smoothing (default: 1)
5. `--mu` - Coulomb pseudopotential (default: 0.1)
6. `--tol` - structure tolerance in angstrom (default: 0.2)

In the directory specified by `-p` the `results/` folder will appear.

Output files are:
1. `direct_s*.csv` file with the a2f lambda, wlog, w2 and Tc calculated directly from lambdas and gammas 
2. `a2f_s*_r*_g*.csv` file with the a2f, lambda, wlog, w2 and Tc calculated by integrating the computed a2f
3. `*.vasp` and `*.cif` file with the crystal structure
4. `plot_s*_r*_g*.pdf` with visualized parameters
5. `plot_article.pdf` article-view plot
6. `result.json` contains all computed parameters of superconducting state

