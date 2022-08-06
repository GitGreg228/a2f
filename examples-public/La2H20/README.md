Input command in the `src/` dir:

```
python main.py -p ..\examples-public\La2H20
```

Output:
~~~
Error: Unable to detect *dyn* files in ..\examples-public\La2H20
Found dyn_elphs in ..\examples-public\La2H20\La2H20.dyn1.elph.1, ..\examples-public\La2H20\La2H20.dyn2.elph.2, ..\examples-public\La2H20\La2H20.dyn3.elph.3, ..\examples-public\La2H20\La2H20.dyn4.elph.4
Found ph.outs in ..\examples-public\La2H20\La2H20.ph.out

The structure under consideration is 2x of Fm-3m-LaH10.
Trying to get weights out of ph.out ...
q = (0.167, 0.167, 0.177) with number of q in the star 8
q = (0.167, -0.500, 0.177) with number of q in the star 4
q = (-0.500, 0.167, 0.177) with number of q in the star 4
q = (-0.500, -0.500, 0.177) with number of q in the star 2
Everything's ok with the weights. They will be used from ph.out file(s).
With the smoothing of 3 THz, direct lambda = 4.004, direct wlog = 736 K, direct w2 = 1112 K.
After calculating alpha2f on resolution r = 1000 and applying gaussian filter with sigma = 5, interp lambda = 4.004, interp wlog = 746 K, interp w2 = 1114 K.
Mu = 0.1. Calculated McMillan Tc = 147.793 (149.571) K, calculated Allen-Dynes Tc = 258.757 (259.697) K
Eliashberg Tc = 276.338+-0.001 K
lambda: 4.004
wlog: 736 K
w2: 1112 K
McM Tc: 147.8 K
AD Tc: 258.8 K
E Tc: 276.3 K
DOS: 0.207 states/spin/Ry/A^3
Sommerfeld gamma: 0.597 mJ/cm^3/K^2
DeltaC/Tc: 0.095 mJ/cm^3/K^2
Delta: 63.2 meV
2Delta/kBTc: 5.31
Hc2: 91.7 T
~~~
