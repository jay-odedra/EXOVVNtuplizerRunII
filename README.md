# EXOVVNtuplizerRunII

Ntuplizer for searches for heavy resonances decaying to dibosons

## installation instructions

Setting up CMSSW (for september reprocessing):

```
cmsrel CMSSW_12_4_3
cd CMSSW_12_4_3/src
cmsenv
git cms-init
```



### getting the ntuplizer code
```
cd $CMSSW_BASE/src
export GITUSER=`git config user.github`
git clone https://github.com/${GITUSER}/EXOVVNtuplizerRunII 
cd EXOVVNtuplizerRunII
git remote add UZHCMS https://github.com/UZHCMS/EXOVVNtuplizerRunII
git fetch UZHCMS
git checkout -b Development_Run3DiElectron UZHCMS/Run3DiElectron
cd $CMSSW_BASE/src
scram b -j 8
```

### running locally
run your config as follows :

```
cd EXOVVNtuplizerRunII/Ntuplizer
cmsRun nanoanalyzercrab_cfg_22_ZeroBias.py
```


### CRAB submission 
do e.g. 
```
crab submit crab.py

```
