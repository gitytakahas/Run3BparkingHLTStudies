# Run3BparkingHLTStudies

Before going through below, you should first setup ROOT env. 

# 1. Setup

```
git clone git@github.com:gitytakahas/Run3BparkingHLTStudies.git
cd Run3BparkingHLTStudies
mkdir root 
cp /afs/cern.ch/user/y/ytakahas/public/forRob/BparkingForRun3/*.root root/
```

# 2. (re-)create efficiency 

```
python compare_eff.py 
```

This will create effmap.root to be later used to make the ROC curve 


# 3. (re-) create ROC curve 

```
python draw_roc.py --pu 50 (e.g. you can change the number here)
```

If you want to have the envelope of the best points (HLT = L1 - 1GeV), then, do, 

```
python draw_roc.py --pu 50 --envelope 
```

The representative points are 
   * pu = 25 (0.9E34)
   * pu = 30 (1.1E34)
   * pu = 36 (1.3E34)
   * pu = 42 (1.5E34)
   * pu = 47 (1.7E34)
   * pu = 56 (2.0E34)
   
   
   

