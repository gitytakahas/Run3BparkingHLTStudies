# Run3BparkingHLTStudies

Unless you want to generate analysis ROOT files by yourself, you can skip Step1 and 2. 

1. Setting up your environemnt https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMHLTRun3RecommendationForPAG

2. Generate EDM file, starting from RAW mc/data. These can be found at 

In case of MC, use /store/user/ytakahas/Trigger/Winter21/*.root 
In case of data, use ```/EphemeralZeroBias*/Run2018D-v1/RAW``` where * stands for [1,8]

3. make flat ntuple, starting from EDM files

> python3 makeRun3Ntup.py 

4. generate ROOT file for the signal efficiency study, using the output in step3.

> python3 eff_hlt.py

5. generate ROOT file for the rate study, using the output in step3.

> python3 rate.py 


6. efficiency plotting using the output in step4.

> python3 compare_eff.py

7. rate plotting using the output in step5.

> python3 compare_rate.py 
