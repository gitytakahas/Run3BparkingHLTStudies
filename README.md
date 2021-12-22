# Run3BparkingHLTStudies

The instructions below have been verified with python3 (CMSSW_12_0_X onwards).

# 1. log-in to lxplus 
```
cmsrel CMSSW_12_1_0
cd CMSSW_12_1_0/src
cmsenv 
```

# 2. Setup package

```
git clone git@github.com:gitytakahas/Run3BparkingHLTStudies.git
cd Run3BparkingHLTStudies
```

# 3. (Can be skipped) derive efficiecny for single-mu after folding in analysis efficiency. 

```
cd single-mu/
python3 draw_roc.py --pr --weight 
```

If you just need L1 times HLT trigger eff. you can do, instead,
```
cd single-mu/
python3 draw_roc.py --pr
```

These numbers have been already embedded into the ```draw_roc.py``` script, which is the one used for drawing the ROC curve.


# 4. (Can be skipped) create efficiency map for di-e after folding in analysis efficiency.

```
python3 compare_eff.py --weight
```

This will create the efficiency map to be later used to make the ROC curve. 
If you don't want to fold in analysis efficiency, you can just remove ```--weight``` option. 
Obtained maps are already stored in eos so this step can be skipped unless you want to create the new one.
If you want to make comparison plots for several basic distributions at the HLT level, just add --plot option.


# 5. (Can be skipped) create rate map for di-electron trigger

```
python3 compare_rate.py
```

This will create the HLT rate map (as a function of number of pileup) to be later used to make the ROC curve. 
Obtained maps are already stored in eos so this step can be skipped unless you want to create the new one.


# 6. (re-) create ROC curve 

```
python3 draw_roc.py --pu 50 --weight (e.g. you can change the number here)
```

If you want to have the envelope of the best points (HLT = L1 - 1GeV), then, do, 

```
python3 draw_roc.py --pu 50 --weight --envelope 
```

The representative points are 
   * pu = 25 (0.9E34)
   * pu = 30 (1.1E34)
   * pu = 36 (1.3E34)
   * pu = 42 (1.5E34)
   * pu = 48 (1.7E34)
   * pu = 56 (2.0E34)

This step will create the ROOT files containing the ROC curve, to be later used in the next step (estimate.py).


# 5. estimate Kee events
```
python3 estimate.py 
```
