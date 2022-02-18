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
python main.py -p path/to/your/folder/with/dyns
```
