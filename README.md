flashgg
=======

Before you start, **please take note** of these warnings and comments:
* **N.B.** Make sure you are on lxplus6 or otherwise using an SLC6 machine. Make sure SCRAM_ARCH is slc6_amd64_gcc491 (or later for 80X).
* **N.B.** The setup script will check out many packages and take a while!
* **N.B.** You can ignore "error: addinfo_cache" lines. 
* **N.B.** This is to set up the latest area in a self-consistent way. If you want a particular flashgg version corresponding to other samples, please see https://twiki.cern.ch/twiki/bin/viewauth/CMS/FLASHggFramework#Instructions_for_users

Get everything you need, starting from a clean area:

 ```
 cmsrel CMSSW_8_0_28
 cd CMSSW_8_0_28/src
 cmsenv
 git cms-init
 cd $CMSSW_BASE/src 
 git clone https://github.com/ArnabPurohit/flashgg flashgg
 cd flashgg
 git checkout forLegacy2016_Working 
 cd ..
 source flashgg/setup.sh
 ```

If everything now looks reasonable, you can build (note: if it does not compile at the first attempt try running "scram b -j 4" once more.):
 ```
 cd $CMSSW_BASE/src
 scram b -j 4
 ```
And a very basic workflow test:
 ```
 cd $CMSSW_BASE/src/flashgg
 cmsRun MicroAOD/test/microAODstd.py processType=sig datasetName=glugluh
 cmsRun Taggers/test/simple_Tag_test.py
 cmsRun Taggers/test/diphotonsDumper_cfg.py
 cmsRun Systematics/test/workspaceStd.py processId=ggh_125
For Validation package
 cmsRun Validation/test/makeTreePhotons_80X.py
 ```

These are just some test examples; the first makes MicroAOD from a MiniAOD file accessed via xrootd, 
the second produces tag objects and screen output from the new MicroAOD file,
and the other two process the MicroAOD file to test ntuple and workspace output.

The setup code will automatically change the initial remote branch's name to upstream to synchronize with the project's old conventions.  
The code will also automatically create an "origin" repo based on its guess as to where your personal flashgg fork is.
Check that this has worked correctly if you have trouble pushing.  (See setup.sh for what it does.)

