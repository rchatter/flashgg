import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

process = cms.Process("tnp")

###################################################################
myoptions = dict()
isMC = False

doShowershape = True

myoptions['HLTProcessName']        = "HLT"
if (doShowershape):
    myoptions['DIPHOTON_COLL']     = "flashggUpdatedIdMVADiPhotons"
else:
    myoptions['DIPHOTON_COLL']     = "flashggDiPhotons"

myoptions['PHOTON_CUTS']           = "(abs(superCluster.eta)<2.5) && ((superCluster.energy*sin(superCluster.position.theta))>20.0)&& !(1.4442<=abs(superCluster.eta)<=1.566)"
myoptions['PHOTON_TAG_CUTS']       = "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && (superCluster.energy*sin(superCluster.position.theta))>40.0"
myoptions['MAXEVENTS']             = cms.untracked.int32(10) # later change to "-1" 
myoptions['useAOD']                = cms.bool(False)
myoptions['OUTPUTEDMFILENAME']     = 'edmFile.root'
myoptions['DEBUG']                 = cms.bool(False)
myoptions['LEADING_PRESELECTION']  = """(abs(leadingPhoton.superCluster.eta) < 2.5 && abs(subLeadingPhoton.superCluster.eta) < 2.5) &&
                                        (leadingPhoton.pt > 30) && (leadingPhoton.hadronicOverEm < 0.08) &&
                                        ((leadingPhoton.full5x5_r9 > 0.5 && leadingPhoton.isEB) || (leadingPhoton.full5x5_r9 > 0.8 && leadingPhoton.isEE)) &&
                                        ((subLeadingPhoton.full5x5_r9 > 0.5 && subLeadingPhoton.isEB) || (subLeadingPhoton.full5x5_r9 > 0.8 && subLeadingPhoton.isEE)) &&
                                        ((leadingPhoton.isEB && 
                                        (leadingPhoton.full5x5_r9>0.85 ||
                                        (leadingPhoton.full5x5_sigmaIetaIeta < 0.015 && leadingPhoton.pfPhoIso03 < 4.0 && leadingPhoton.trkSumPtHollowConeDR03 < 6.0 ))) ||
                                        (leadingPhoton.isEE && 
                                        (leadingPhoton.full5x5_r9>0.9 ||
                                        (leadingPhoton.full5x5_sigmaIetaIeta < 0.035 && leadingPhoton.pfPhoIso03 < 4.0 && leadingPhoton.trkSumPtHollowConeDR03 < 6.0 )))) &&
                                        (leadingPhoton.pt > 14 && leadingPhoton.hadTowOverEm()<0.15 &&
                                        (leadingPhoton.full5x5_r9()>0.8 || leadingPhoton.chargedHadronIso()<20 || leadingPhoton.chargedHadronIso()<0.3*leadingPhoton.pt())) 
""" 

myoptions['SUBLEADING_PRESELECTION'] = """(abs(leadingPhoton.superCluster.eta) < 2.5 && abs(subLeadingPhoton.superCluster.eta) < 2.5) &&
                                          (subLeadingPhoton.pt > 20) && (subLeadingPhoton.hadronicOverEm < 0.08) &&
                                          ((leadingPhoton.full5x5_r9 > 0.5 && leadingPhoton.isEB) || (leadingPhoton.full5x5_r9 > 0.8 && leadingPhoton.isEE)) &&
                                          ((subLeadingPhoton.full5x5_r9 > 0.5 && subLeadingPhoton.isEB) || (subLeadingPhoton.full5x5_r9 > 0.8 && subLeadingPhoton.isEE)) &&
                                          (( subLeadingPhoton.isEB && 
                                          (subLeadingPhoton.full5x5_r9>0.85 ||
                                          (subLeadingPhoton.full5x5_sigmaIetaIeta < 0.015 && subLeadingPhoton.pfPhoIso03 < 4.0 && subLeadingPhoton.trkSumPtHollowConeDR03 < 6.0 ))) ||
                                          (subLeadingPhoton.isEE && 
                                          (subLeadingPhoton.full5x5_r9>0.9 ||
                                          (subLeadingPhoton.full5x5_sigmaIetaIeta < 0.035 && subLeadingPhoton.pfPhoIso03 < 6.0 && subLeadingPhoton.trkSumPtHollowConeDR03 < 6.0 )))) &&
                                          (subLeadingPhoton.pt > 14 && subLeadingPhoton.hadTowOverEm()<0.15 &&
                                          (subLeadingPhoton.full5x5_r9()>0.8 || subLeadingPhoton.chargedHadronIso()<20 || subLeadingPhoton.chargedHadronIso()<0.3*subLeadingPhoton.pt()))
"""


from flashgg.Validation.treeMakerOptionsPhotons_cfi import *

if (isMC):
    myoptions['INPUT_FILE_NAME']       = ("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/170212_190415/0000/myMicroAODOutputFile_101.root")

    myoptions['OUTPUT_FILE_NAME']      = "TnPTree_mc.root"
    myoptions['TnPPATHS']              = cms.vstring("HLT_Ele27_WPTight_Gsf_v*") #HLT_Ele23_WPLoose_Gsf_v #HLT_Ele27_WPTight_Gsf_v
    myoptions['TnPHLTTagFilters']      = cms.vstring("hltEle27WPTightGsfTrackIsoFilter") #hltEle23WPLooseGsfTrackIsoFilter  #hltEle27WPTightGsfTrackIsoFilter 
    myoptions['TnPHLTProbeFilters']    = cms.vstring()
    myoptions['HLTFILTERTOMEASURE']    = cms.vstring("")
    myoptions['GLOBALTAG']             = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
    myoptions['EVENTSToPROCESS']       = cms.untracked.VEventRange()
else:
    #myoptions['INPUT_FILE_NAME']       = ("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/LegacyReReco-18Apr2017-2_6_0/2_6_0/SingleElectron/LegacyReReco-18Apr2017-2_6_0-2_6_0-v0-Run2016C-18Apr2017-v1/170526_151510/0000/myMicroAODOutputFile_1.root")
    myoptions['INPUT_FILE_NAME']       = ("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/LegacyReReco-07Aug2017-2_6_1/2_6_1/SingleElectron/LegacyReReco-07Aug2017-2_6_1-2_6_1-v0-Run2016C-07Aug17-v1/171008_210948/0000/myMicroAODOutputFile_10.root")
    myoptions['OUTPUT_FILE_NAME']      = "TnPTree_data.root"
    myoptions['TnPPATHS']              = cms.vstring("HLT_Ele27_WPTight_Gsf_v*")
    myoptions['TnPHLTTagFilters']      = cms.vstring("hltEle27WPTightGsfTrackIsoFilter")
    myoptions['TnPHLTProbeFilters']    = cms.vstring()
    myoptions['HLTFILTERTOMEASURE']    = cms.vstring("")
    myoptions['GLOBALTAG']             = '80X_dataRun2_2016SeptRepro_v3'
    myoptions['EVENTSToPROCESS']       = cms.untracked.VEventRange()

###################################################################

setModules(process, myoptions)
from flashgg.Validation.treeContentPhotons_cfi import *

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = myoptions['GLOBALTAG']

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(myoptions['INPUT_FILE_NAME']),
                            eventsToProcess = myoptions['EVENTSToPROCESS']
                            )

process.maxEvents = cms.untracked.PSet( input = myoptions['MAXEVENTS'])

###################################################################
## ID
###################################################################

from flashgg.Validation.photonIDModules_cfi import *
setIDs(process, myoptions)

###################################################################
## SEQUENCES
###################################################################

process.egmPhotonIDs.physicsObjectSrc = cms.InputTag("photonFromDiPhotons")

process.pho_sequence = cms.Sequence(
    process.photonFromDiPhotons +
    process.goodPhotonTags +
    process.goodPhotonProbes +
    process.goodPhotonProbesPreselection +
    process.goodPhotonProbesIDMVA +
    process.goodPhotonTagsIDMVA +
    process.goodPhotonsTagHLT +
    process.goodPhotonsProbeHLT #+ #26.07.2017
    # process.goodPhotonProbesL1
    )

###################################################################
## TnP PAIRS
###################################################################

process.allTagsAndProbes = cms.Sequence()
process.allTagsAndProbes *= process.tagTightRECO


##########################################################################
## TREE MAKER OPTIONS
##########################################################################
if (not isMC):
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )
    
process.PhotonToRECO = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                      mcTruthCommonStuff, CommonStuffForPhotonProbe,
                                      tagProbePairs = cms.InputTag("tagTightRECO"),
                                      arbitration   = cms.string("None"),
                                      flags         = cms.PSet(passingPresel  = cms.InputTag("goodPhotonProbesPreselection"),
                                                               passingIDMVA   = cms.InputTag("goodPhotonProbesIDMVA"),
                                                               ),                                               
                                      allProbes     = cms.InputTag("goodPhotonsProbeHLT"),
                                      )

if (isMC):
    # find out which pair collection is used in PhotonToRECO
    tagProbeCollName = process.PhotonToRECO.tagProbePairs.getModuleLabel()

    print "tagProbeCollName", tagProbeCollName # "tagTightRECO" is the output
    tagProbeCollModule = getattr(process,tagProbeCollName)

    # get collections from which tag/probe pairs were build
    # (assuming they are build in this process here)
    import re
    tagColl, probeColl = re.split('\s+', tagProbeCollModule.decay.value().strip())

    #----------
    # define modules to calculate the matching with generator level ELECTRONs
    #----------

    # adapted from EgammaAnalysis/TnPTreeProducer/python/egmTreesSetup_cff.py
    cut_gen_standard = 'abs(pdgId) == 11 && pt > 3 && abs(eta) < 2.7 && isPromptFinalState'
    cut_gen_flashgg  = 'abs(pdgId) == 11 && pt > 3 && abs(eta) < 2.7 && ( isPromptFinalState || status == 23)'
    cut_gen_tau      = 'abs(pdgId) == 11 && pt > 3 && abs(eta) < 2.7 && ( isPromptFinalState || isDirectPromptTauDecayProductFinalState) '

    process.genEle   = cms.EDFilter( "GenParticleSelector",
                                     src = cms.InputTag('flashggPrunedGenParticles'),
                                     cut = cms.string(cut_gen_standard),
                                     )

    for matchedColl, recoColl in (
        ('genTagPho', tagColl),
        ('genProbePho', probeColl),
        ):
        setattr(process, matchedColl, cms.EDProducer("MCMatcher",
                                                     src      = cms.InputTag(recoColl),
                                                     matched  = cms.InputTag('genEle'), # match reco photons to generated electrons
                                                     mcStatus = cms.vint32(),
                                                     mcPdgId  = cms.vint32(),
                                                     checkCharge = cms.bool(False),
                                                     maxDeltaR   = cms.double(0.20),   # Minimum deltaR for the match
                                                     maxDPtRel   = cms.double(50.0),    # Minimum deltaPt/Pt for the match
                                                     resolveAmbiguities    = cms.bool(False), # Forbid two RECO objects to match to the same GEN objec
                                                     resolveByMatchQuality = cms.bool(True),  # False = just match input in order; True = pick lowest deltaR pair first
                                                     ))
    
    # build the MC sequence
    process.mc_sequence = cms.Sequence(process.genEle + 
                                       process.genTagPho + 
                                       process.genProbePho
                                       )

    process.PhotonToRECO.tagMatches    = cms.InputTag("genTagPho")
    process.PhotonToRECO.probeMatches  = cms.InputTag("genProbePho")

    process.PhotonToRECO.eventWeight   = cms.InputTag("generator")
    # process.PhotonToRECO.eventWeight   = cms.double(1.0) #just to keep 1.0 nothing change in "totWeight" if we put 2.0
    process.PhotonToRECO.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

else:
    # empty MC sequence
    process.mc_sequence = cms.Sequence()


    

process.tree_sequence = cms.Sequence(process.PhotonToRECO)

##########################################################################
## PATHS
##########################################################################

process.out = cms.OutputModule("PoolOutputModule", 
                               fileName = cms.untracked.string(myoptions['OUTPUTEDMFILENAME']),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
                               )
process.outpath = cms.EndPath(process.out)
if (not myoptions['DEBUG']):
    process.outpath.remove(process.out)

##########################################################################################
###### MICROAOD STUFF
##########################################################################################

if (doShowershape):
    process.load("flashgg/Taggers/flashggUpdatedIdMVADiPhotons_cfi")
    process.flashggUpdatedIdMVADiPhotons.DiPhotonTag = cms.InputTag('flashggDiPhotons') # To include shower shape corrections
    process.flashggUpdatedIdMVADiPhotons.doIsoCorrection = cms.bool(False)
    #process.flashggUpdatedIdMVADiPhotons.isoCorrectionFile = cms.FileInPath("flashgg/Taggers/data/pho_iso_corrections_hybrid_LegacyReReco2016_v1.root")
    process.flashggUpdatedIdMVADiPhotons.doNon5x5transformation = cms.bool(False)
    process.flashggUpdatedIdMVADiPhotons.do5x5correction = cms.bool(False)
    #process.flashggUpdatedIdMVADiPhotons.correctionFile = cms.FileInPath("flashgg/MicroAOD/data/transformationforProbe_legacy2016.root")

    process.load("flashgg/Taggers/flashggDiPhotonMVA_cfi")
    # process.flashggDiPhotonMVA.DiPhotonTag = cms.InputTag('flashggDiPhotons')
    process.flashggDiPhotonMVA.DiPhotonTag = cms.InputTag('flashggUpdatedIdMVADiPhotons')
else:
    process.load("flashgg/Taggers/flashggDiPhotonMVA_cfi")
    process.flashggDiPhotonMVA.DiPhotonTag = cms.InputTag('flashggDiPhotons')

# DEBUG: uncomment to disable HLT criterion
# process.hltFilter= cms.EDFilter("HLTBool",
#                                 result = cms.bool(True))


if (doShowershape):

    if (isMC):
        process.p = cms.Path(
            process.flashggUpdatedIdMVADiPhotons +
            process.flashggDiPhotonMVA +
            #process.sampleInfo + # gkole 2nd Jan
            process.hltFilter +
            process.pho_sequence + 
            process.allTagsAndProbes +
            process.pileupReweightingProducer +
            process.mc_sequence + 
            process.tree_sequence
            )
    else:
        process.p = cms.Path(
            process.flashggUpdatedIdMVADiPhotons +
            process.flashggDiPhotonMVA +
            #process.sampleInfo +
            process.hltFilter +
            process.pho_sequence + 
            process.allTagsAndProbes +
            process.pileupReweightingProducer +
            process.mc_sequence +
            process.tree_sequence
            )

else:
    if (isMC):
        process.p = cms.Path(
            process.flashggDiPhotonMVA +
            # process.sampleInfo + # gkole 2nd Jan
            process.hltFilter +
            process.pho_sequence + 
            process.allTagsAndProbes +
            process.pileupReweightingProducer +
            process.mc_sequence + 
            process.tree_sequence
            )
    else:
        process.p = cms.Path(
            process.flashggDiPhotonMVA +
            #process.sampleInfo +
            process.hltFilter +
            process.pho_sequence + 
            process.allTagsAndProbes +
            process.pileupReweightingProducer +
            process.mc_sequence +
            process.tree_sequence
            )
        

process.TFileService = cms.Service("TFileService", 
                                   #fileName = cms.string(myoptions['OUTPUT_FILE_NAME']),
                                   fileName = cms.string("TnP.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

from flashgg.MetaData.JobConfig import JobConfig

customize = JobConfig(crossSections=["$CMSSW_BASE/src/flashgg/MetaData/data/cross_sections.json"])
customize.setDefault("maxEvents", 100)
customize.setDefault("targetLumi", 1)
customize(process)

# get filters summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
