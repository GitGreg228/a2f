Input command in the `src/` dir:

```
python main.py -p ../examples-public/La2H20 -r 1000 -s 6 -g 5
```

Output:
```
Found dyn_elphs in ../examples-public/La2H20\La2H20.dyn1.elph.1, ../examples-public/La2H20\La2H20.dyn2.elph.2, ../examples-public/La2H20\La2H20.dyn3.elph.3, ../examples-public/La2H20\La2H20.dyn4.elph.4
Found ph.outs in ../examples-public/La2H20\La2H20.ph.out

The structure under consideration is Fm-3m-LaH10.
q = (0.167, 0.167, 0.177) with number of q in the star 8
q = (0.167, -0.500, 0.177) with number of q in the star 4
q = (-0.500, 0.167, 0.177) with number of q in the star 4
q = (-0.500, -0.500, 0.177) with number of q in the star 2
With the smoothing of 6 THz, direct lambda = 3.184, direct wlog = 902 K, direct w2 = 1220 K.
After calculating alpha2f on resolution r = 1000 and applying gaussian filter with sigma = 5, interp lambda = 3.186, interp wlog = 924 K, interp w2 = 1229 K.
Mu = 0.1. Calculated McMillan Tc = 166.468 (170.625) K, calculated Allen-Dynes Tc = 245.383 (245.383) K
Eliashberg Tc = 265.305+-0.001 K
```
