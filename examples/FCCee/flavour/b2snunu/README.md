# How to run
## Existing method for building since this branch has already been updated to be in line with pre-edm4hep1, execute in head directory for FCCAnalyses
```bash
source setup.sh
fccanalysis build -j 16
```

## Run stage1 of analysis
In general avoid using ncpus>1 as I suspect there may be a race condition in the vertexing, and the speed doesn't scale very much anyway
```bash
cd examples/FCCee/flavour/b2snunu/
# test
fccanalysis run analysis_stage1.py --test --decay Lb2LNuNu --nevents 1000
# full run, single job (using fraction to reduce portion of input used)
fccanalysis run analysis_stage1.py --decay Lb2LNuNu --fraction 0.0001 --output-dir /eos/lhcb/user/a/aiwieder/FCC/b2snunu/output/stage1 |& tee stage1.log
# run on condor, the quickest a job on a single qqbar file can be is ~40 minutes assuming a standard rate of 40 events per second -> chunks 1000, but some speed up to be had by running over more than one file so can play around with fewer chunks
# below is just an initial guess at how we might want to run it, can use --fraction just to run on a small portion again
fccanalysis run analysis_stage1.py --decay Bd2KstNuNu --output-dir /eos/lhcb/user/a/aiwieder/FCC/b2snunu/output/stage1/ --condor --chunks 500 --ncpus 1 |& tee /eos/lhcb/user/a/aiwieder/FCC/b2snunu/output/logs/stage1_submit.log
# condor logs get written to `FCCAnalyses/BatchOutputs/`
```

# Notes
## original building with pre-edm4hep1 so we can use old MC
```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
git clone --branch pre-edm4hep1 git@github.com:HEP-FCC/FCCAnalyses.git
cd FCCAnalyses
source ./setup.sh
fccanalysis build -j 8
```