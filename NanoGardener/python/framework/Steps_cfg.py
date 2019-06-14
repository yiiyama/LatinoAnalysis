#!/usr/bin/env python

## --------------------------------- Some predefined sequence Chains -----------------------------------------

#PU2017Seq=['puW2017B','puW2017C','puW2017D','puW2017E','puW2017F','puW2017perEra']   
PU2017Seq=['puW']

# -------------------------------------------- HERE WE GO ----------------------------------------------------

Steps = {

# ------------------------------------------------ CHAINS ----------------------------------------------------

  'puW2017Test': {
                  'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'subTargets' : ['trigMC'] + PU2017Seq ,
                 }, 

  'MCnofilter' : {
                  'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'subTargets' : ['leptonMaker'],
                 },

## ------- MC:

  'MCl1loose2016': {
                  'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'selection'  : '"(nElectron>0 && Electron_pt[0]>10) || (nMuon>0 && Muon_pt[0]>10)"' , 
                  'subTargets' : ['baseW', 'leptonMaker','lepSel', 'puW2016', 'l2Kin', 'l3Kin', 'l4Kin', 'btagPerJet2016', 'btagPerEvent'],
                },

  'MCl1loose2017': {
                  'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'selection'  : '"((nElectron+nMuon)>0)"' ,
                  'subTargets' : ['leptonMaker','lepSel', 'puW2017', 'l2Kin', 'l3Kin', 'l4Kin', 'btagPerJet2017', 'btagPerEvent'],
                },

  'MCl1loose2017v2': {
                  'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'selection'  : '"((nElectron+nMuon)>0)"' ,
                  'subTargets' : ['leptonMaker','lepSel', 'puW2017', 'l2Kin', 'l3Kin', 'l4Kin', 'btagPerJet2017', 'btagPerEvent','PrefCorr2017'],
                },

  'MCCorr2017' : {
                 'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'subTargets' : ['baseW','PrefCorr2017','jetSel','CleanJetCut',
                                  'PromptParticlesGenVars','GenVar','GenLeptonMatch', 'HiggsGenVars', 'TopGenVars', 'wwNLL',
                                  'rochesterMC','trigMC','LeptonSF','puW','l2Kin', 'l3Kin', 'l4Kin','formulasMC'],
                    },


  'MCCorr2017test' : { 
                 'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'subTargets' : ['baseW','PrefCorr2017','jetSel',
                                  'PromptParticlesGenVars','GenVar','GenLeptonMatch', 'HiggsGenVars', 'TopGenVars', 'wwNLL', 
                                  'trigMC','LeptonSF','puW','l2Kin', 'l3Kin', 'l4Kin','formulasMC'], 
                   'onlySample' : [
                                  #### DY
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1','DYJetsToLL_M-5to50-LO',
                                  ####
                                # 'WJetsToLNu','WJetsToLNu_HT100_200','WJetsToLNu_HT200_400','WJetsToLNu_HT400_600','WJetsToLNu_HT600_800',
                                # 'WJetsToLNu_HT800_1200','WJetsToLNu_HT1200_2500','WJetsToLNu_HT2500_inf',
                                  ####
                                  'QCD_Pt-15to20_EMEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched','QCD_Pt-50to80_EMEnriched_ext1',
                                  'QCD_Pt-20toInf_MuEnrichedPt15','QCD_Pt-30toInf_DoubleEMEnriched','QCD_Pt-15to20_MuEnrichedPt5',
                                  ####
                                  'QCD_Pt_15to20_bcToE','QCD_Pt_20to30_bcToE','QCD_Pt_30to80_bcToE','QCD_Pt_80to170_bcToE',
                                  'QCD_Pt_170to250_bcToE','QCD_Pt_250toInf_bcToE',
                                  ####
                                  'TT','TTJets','TTTo2L2Nu', 
                                  ]  
                },


  'MCformulas': {
                  'isChain'    : True  ,
                  'do4MC'      : True ,
                  'do4Data'    : False  ,           
                  'subTargets' : ['baseW'],#,'trigMC','formulasMC'],
                },             

  'MCWeights2017' : {
                  'isChain'    : True  ,
                  'do4MC'      : True ,
                  'do4Data'    : False  ,
                  'subTargets' : ['PrefCorr2017','jetSel','PromptParticlesGenVars','GenVar','GenLeptonMatch', 'HiggsGenVars', 'TopGenVars', 'wwNLL', 'trigMC','LeptonSF','formulasMC'],
                  #'subTargets' : ['baseW','trigMC','LeptonSF','formulasMC'],
                },             




## ------- WgStar MC:

  'MCWgStar2017' : { 
                     'isChain'    : True  ,
                     'do4MC'      : True  ,
                     'do4Data'    : False ,
                     'selection'  : '"((nElectron+nMuon)>1)"' ,
                     'subTargets' : ['leptonMaker','WgSSel','puW2017', 'l2Kin', 'l3Kin', 'l4Kin', 'btagPerJet2017', 'btagPerEvent',
                                     'PromptParticlesGenVars','GenVar','GenLeptonMatch', 'HiggsGenVars', 'TopGenVars', 'wwNLL', 'trigMC','LeptonSF'],
                     'onlySample' : [
                                   'Wg500','Wg_AMCNLOFXFX','WZTo3LNu','Wg_MADGRAPHMLM',
                                   #'Wg500','Wg_AMCNLOFXFX','WZTo3LNu','WgStarLNuEE','WgStarLNuMuMu','Wg_MADGRAPHMLM',
                                   'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3',
                                   'DYJetsToLL_M-5to50-LO','DYJetsToLL_M-50-LO-ext1',
                                   'WZTo2L2Q','WZTo3LNu_mllmin01_ext1','WZTo3LNu',
                                 ]
                   },

## ------- DATA:
    
  'DATAl1loose2016': {
                  'isChain'    : True  ,
                  'do4MC'      : False ,
                  'do4Data'    : True  ,
                  'selection'  : '"(nElectron>0 && Electron_pt[0]>10) || (nMuon>0 && Muon_pt[0]>10)"' , 
                  'subTargets' : ['leptonMaker','lepSel', 'l2Kin', 'l3Kin', 'l4Kin'],
                },

  'DATAl1loose2017': {
                  'isChain'    : True  ,
                  'do4MC'      : False ,
                  'do4Data'    : True  ,
                  'selection'  : '"((nElectron+nMuon)>0)"' ,
                  'subTargets' : ['leptonMaker','lepSel', 'l2Kin', 'l3Kin', 'l4Kin','trigData','formulasDATA'],
                }, 

  'DATAl1loose2017v2': {
                  'isChain'    : True  ,
                  'do4MC'      : False ,
                  'do4Data'    : True  ,
                  'selection'  : '"((nElectron+nMuon)>0)"' ,
                  'subTargets' : ['leptonMaker','lepSel', 'l2Kin', 'l3Kin', 'l4Kin','trigData','formulasDATA'],
                },

  'DATAFakeCorr2017': {
                  'isChain'    : True  ,
                  'do4MC'      : False ,
                  'do4Data'    : True  ,
                  'subTargets' : ['trigData','rochesterDATA','jetSel','l2Kin','formulasDATA'],
                },

  'DATACorr2017' : {
                  'isChain'    : True  ,
                  'do4MC'      : False ,
                  'do4Data'    : True  ,
                  'subTargets' : ['rochesterDATA','jetSel','l2Kin', 'l3Kin', 'l4Kin','formulasDATA'],
                },

  'DATACorr2017step2' : {
                  'isChain'    : True  ,
                  'do4MC'      : False ,
                  'do4Data'    : True  ,
                  'subTargets' : ['CleanJetCut','l2Kin','l3Kin', 'l4Kin'],
                },

## ------- WgStar DATA:

    'DATAWgStar2017' : { 
                  'isChain'    : True  ,
                  'do4MC'      : False ,
                  'do4Data'    : True  ,
                  'selection'  : '"((nElectron+nMuon)>1)"' ,
                  'subTargets' : ['leptonMaker','WgSSel', 'l2Kin', 'l3Kin', 'l4Kin','trigData','formulasDATA'],
                   },

#  Merged back to DATAl1loose2017 (was a tmp fix)    
#  'DATAformulas' : {
#                 'isChain'    : True  ,
#                 'do4MC'      : False ,
#                 'do4Data'    : True  ,           
#                 'subTargets' : ['trigData','formulasDATA'],
#               },

# ------------------------------------------------ MODULES ---------------------------------------------------

## ------- MODULES: MC Kinematic
  
  'PromptParticlesGenVars' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.PromptParticlesGenVarsProducer' ,
                  'declare'    : 'PromptParticlesGenVars = lambda : PromptParticlesGenVarsProducer()',
                  'module'     : 'PromptParticlesGenVars()',
                  } , 


  'GenVar'       : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.GenVarProducer' ,
                  'declare'    : 'GenVar = lambda : GenVarProducer()',
                  'module'     : 'GenVar()' ,
                   },

  'GenLeptonMatch' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.GenLeptonMatchProducer' ,
                  'declare'    : 'GenLeptonMatch = lambda : GenLeptonMatchProducer()',
                  'module'     : 'GenLeptonMatch()' ,
                   },

   'HiggsGenVars' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.HiggsGenVarsProducer' ,
                  'declare'    : 'HiggsGenVars = lambda : HiggsGenVarsProducer()',
                  'module'     : 'HiggsGenVars()',
                  } ,                 

   'TopGenVars' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.TopGenVarsProducer' ,
                  'declare'    : 'TopGenVars = lambda : TopGenVarsProducer()',
                  'module'     : 'TopGenVars()',
                  'onlySample' : ['TTTo2L2Nu', 'TTToSemileptonic']
                  } ,

    'wwNLL' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.wwNLLcorrectionWeightProducer' ,
                  'declare'    : 'wwNLL = lambda : wwNLLcorrectionWeightProducer()',
                  'module'     : 'wwNLL()',
                  'onlySample' : ['WW-LO', 'WWTo2L2Nu', 'WWTo2L2Nu_CP5Up', 'WWTo2L2Nu_CP5Down']
                  } ,

## ------- MODULES: Object Handling

  'lepMergerHWW' : { 
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger' ,
                  'declare'    : 'lepMergerHWW = lambda : collectionMerger( input  = ["Electron","Muon"], output = "Lepton", reverse = True)' ,
                  'module'     : 'lepMergerHWW()' ,
               },  

  'leptonMaker': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.LeptonMaker' ,
                  'declare'    : 'leptonMaker = lambda : LeptonMaker()' ,
                  'module'     : 'leptonMaker()' ,
               }, 

   'lepSel': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.LeptonSel' ,
                  'declare'    : 'leptonSel = lambda : LeptonSel("RPLME_CMSSW", "Loose", 1)' ,
                  'module'     : 'leptonSel()' ,
               },

   'WgSSel' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.LeptonSel' ,
                  'declare'    : 'leptonSel = lambda : LeptonSel("RPLME_CMSSW", "WgStar", 2)' ,
                  'module'     : 'leptonSel()' ,
               },             

   'jetSel'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.JetSel' ,
                  # jetid=2,pujetid='loose',minpt=15.0,maxeta=4.7,jetColl="CleanJet"
                  'declare'    : 'jetSel = lambda : JetSel(2,"loose",15.0,4.7,"CleanJet")' ,
                  'module'     : 'jetSel()' ,
               }, 

   'CleanJetCut' : {
                 'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.CopyCleanJet',
                  'declare'    : 'cleanJetCut = lambda : CopyCleanJet(newcollectionname="CleanJetCut", cuts=["eta>2.65","eta<3.139"])',
                  'module'     : 'cleanJetCut()',
               }, 


## ------- MODULES: Trigger

  'PrefCorr2017' : { 
                 'isChain'    : False ,
                 'do4MC'      : True ,
                 'do4Data'    : False  ,
                 'import'     : 'LatinoAnalysis.NanoGardener.modules.PrefireCorr' ,
                 'declare'    : 'prefCorr2017 = lambda : PrefCorr(jetroot="L1prefiring_jetpt_2017BtoF.root", jetmapname="L1prefiring_jetpt_2017BtoF", photonroot="L1prefiring_photonpt_2017BtoF.root", photonmapname="L1prefiring_photonpt_2017BtoF", variation=0, UseEMpT=0)',
                 'module'     : 'prefCorr2017()',
               },

  'trigData' : { 'isChain'    : False ,
                 'do4MC'      : False ,
                 'do4Data'    : True  ,
                 'import'     : 'LatinoAnalysis.NanoGardener.modules.TrigMaker' ,
                 'declare'    : 'trigData = lambda : TrigMaker("RPLME_CMSSW",True)',
                 'module'     : 'trigData()',
               },

 
  'trigMC'   : { 'isChain'    : False ,
                 'do4MC'      : True  ,
                 'do4Data'    : False ,
                 'import'     : 'LatinoAnalysis.NanoGardener.modules.TrigMaker' ,
                 'declare'    : 'trigMC = lambda : TrigMaker("RPLME_CMSSW",False,False)',
                 'module'     : 'trigMC()',
               },

  'trigMCKeepRun' : { 'isChain'    : False ,
                 'do4MC'      : True  ,
                 'do4Data'    : False ,
                 'import'     : 'LatinoAnalysis.NanoGardener.modules.TrigMaker' ,
                 'declare'    : 'trigMCKR = lambda : TrigMaker("RPLME_CMSSW",False,True)',
                 'module'     : 'trigMCKR()',
               },

## ------- MODULES: MC Weights

  'baseW'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.Grafter' ,
                  'module'     : 'Grafter(["baseW/F=RPLME_baseW","Xsec/F=RPLME_XSection"])',
               },  

  'wwNLL'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.wwNLLcorrectionWeightProducer' ,
                  'module)'     : 'wwNLLcorrectionWeightProducer()',  
               },  

  'btagPerJet2016': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer' ,
                  'module'     : 'btagSFProducer(era="2016", algo="cmva")',
                 },

  'btagPerJet2017': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer' ,
                  'declare'    : 'btagSFProducer2017 = lambda : btagSFProducer(era="2017", algo="deepcsv")',
                  'module'     : 'btagSFProducer2017()',
                 },               

  'btagPerEvent': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.BTagEventWeightProducer' ,
                  'declare'    : '',
                  'module'     : 'BTagEventWeightProducer()',
        
                },


  'LeptonSF' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.LeptonSFMaker' ,
                  'declare'    : 'LeptonSF = lambda : LeptonSFMaker("RPLME_CMSSW")',
                  'module'     : 'LeptonSF()',
                   },

## ------- Pile-Up weights

  'puW'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.runDependentPuW' ,
                  'declare'    : 'puWeight = lambda : runDependentPuW("RPLME_CMSSW")',
                  'module'     : 'puWeight()', 
             } , 

  'puW2016': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer' ,
                  'declare'    : 'pufile_mc2016="%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/pileup/pileup_profile_Summer16.root" % os.environ["CMSSW_BASE"]; pufile_data2016="%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/pileup/PileupData_GoldenJSON_Full2016.root" % os.environ["CMSSW_BASE"]',
                  'module'     : 'puWeightProducer(pufile_mc2016,pufile_data2016,"pu_mc","pileup",verbose=False)',

                },
              
  'puW2017': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer' ,
                  'declare'    : 'pufile_data2017="%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/pileup/pileup_Cert_294927-306462_13TeV_PromptReco_Collisions17_withVar.root" % os.environ["CMSSW_BASE"]',
                  'module'     : 'puWeightProducer("auto",pufile_data2017,"pu_mc","pileup",verbose=False)',
  },

  'puW2017B': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer' ,
                  'declare'    : 'pufile_data2017B="%s/src/LatinoAnalysis/NanoGardener/python/data/PUweights/2017/PU_2017B.root" % os.environ["CMSSW_BASE"]',
                  'module'     : 'puWeightProducer("auto",pufile_data2017B,"pu_mc","pileup",verbose=False,name="puW2017B",doSysVar=False)',
  },


  'puW2017C': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer' ,
                  'declare'    : 'pufile_data2017C="%s/src/LatinoAnalysis/NanoGardener/python/data/PUweights/2017/PU_2017C.root" % os.environ["CMSSW_BASE"]',
                  'module'     : 'puWeightProducer("auto",pufile_data2017C,"pu_mc","pileup",verbose=False,name="puW2017C",doSysVar=False)',
  },


  'puW2017D': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer' ,
                  'declare'    : 'pufile_data2017D="%s/src/LatinoAnalysis/NanoGardener/python/data/PUweights/2017/PU_2017D.root" % os.environ["CMSSW_BASE"]',
                  'module'     : 'puWeightProducer("auto",pufile_data2017D,"pu_mc","pileup",verbose=False,name="puW2017D",doSysVar=False)',
  },


  'puW2017E': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer' ,
                  'declare'    : 'pufile_data2017E="%s/src/LatinoAnalysis/NanoGardener/python/data/PUweights/2017/PU_2017E.root" % os.environ["CMSSW_BASE"]',
                  'module'     : 'puWeightProducer("auto",pufile_data2017E,"pu_mc","pileup",verbose=False,name="puW2017E",doSysVar=False)',
  },

  'puW2017F': {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer' ,
                  'declare'    : 'pufile_data2017F="%s/src/LatinoAnalysis/NanoGardener/python/data/PUweights/2017/PU_2017F.root" % os.environ["CMSSW_BASE"]',
                  'module'     : 'puWeightProducer("auto",pufile_data2017F,"pu_mc","pileup",verbose=False,name="puW2017F",doSysVar=False)',
  },

  'puW2017perEra' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.GenericFormulaAdder' ,
                  'declare'    : '',
                  'module'     : 'GenericFormulaAdder(\'data/formulasToAdd_puW2017.py\')' ,
                 },


## ------- MODULES: Fakes

  'fakeW'  : {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'subTargets' : ['fakeWstep','formulasFAKE'],
                   },

  'fakeWstep'   : {
                  'isChain'    : False ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.LeptonFakeWMaker',
                  'declare'    : '',
                  'module'     : 'LeptonFakeWMaker("RPLME_CMSSW")',
              },

    'fakeSkimDATA': {
        'isChain': False,
        'do4MC': False,
        'do4Data': True,
        'import': 'LatinoAnalysis.NanoGardener.modules.FakeSkimMaker',
        'module': 'FakeSkimMaker("RPLME_CMSSW", ["OppositeSignDielectronSkim", "DimuonElectronSkim", "DimuonPhotonSkim", "SameSignMuonElectronSkim"], False, False)'
    },

    'fakeSkimMC': {
        'isChain': False,
        'do4MC': True,
        'do4Data': False,
        'import': 'LatinoAnalysis.NanoGardener.modules.FakeSkimMaker',
        'module': 'FakeSkimMaker("RPLME_CMSSW", ["OppositeSignDielectronSkim", "DimuonElectronSkim", "DimuonPhotonSkim", "SameSignMuonElectronSkim"], True, False)'
    },

    'fakeSkimJetElectronDATA': {
        'isChain': False,
        'do4MC': False,
        'do4Data': True,
        'import': 'LatinoAnalysis.NanoGardener.modules.FakeSkimMaker',
        'module': 'FakeSkimMaker("RPLME_CMSSW", ["JetElectronSkim", "SameSignDielectronSkim"], False, False)',
        'onlySample': [
            "SingleElectron_Run2017B-31Mar2018-v1",
            "SingleElectron_Run2017C-31Mar2018-v1",
            "SingleElectron_Run2017D-31Mar2018-v1",
            "SingleElectron_Run2017E-31Mar2018-v1",
            "SingleElectron_Run2017F-31Mar2018-v1",
            "DoubleEG_Run2017B-31Mar2018-v1",
            "DoubleEG_Run2017C-31Mar2018-v1",
            "DoubleEG_Run2017D-31Mar2018-v1",
            "DoubleEG_Run2017E-31Mar2018-v1",
            "DoubleEG_Run2017F-31Mar2018-v1",
        ]
    },

    'fakeSkimJetElectronMC': {
        'isChain': False,
        'do4MC': True,
        'do4Data': False,
        'import': 'LatinoAnalysis.NanoGardener.modules.FakeSkimMaker',
        'module': 'FakeSkimMaker("RPLME_CMSSW", ["JetElectronSkim", "SameSignDielectronSkim", "ElectronSkim"], True, False)'
    },

## ------- MODULES: Rochester corrections

  'rochesterMC'   : {
                  'isChain'    : False ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.rochester_corrections',
                  'declare'    : 'rochesterMC = lambda : rochester_corr(False)',
                  'module'     : 'rochesterMC()',
              },

  'rochesterDATA'   : {
                  'isChain'    : False ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.rochester_corrections',
                  'declare'    : 'rochesterDATA = lambda : rochester_corr(True)',
                  'module'     : 'rochesterDATA()',
              },

## ------- MODULES: Kinematic

  'l2Kin'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.l2KinProducer' ,
                  'declare'    : '',
                  'module'     : 'l2KinProducer()' ,
               },  

  'l3Kin'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.l3KinProducer' ,
                  'declare'    : '',
                  'module'     : 'l3KinProducer()' ,
               },  

  'l4Kin'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.l4KinProducer' ,
                  'declare'    : '',
                  'module'     : 'l4KinProducer()' ,
               },  

## ------- MODULES: Adding Formulas

  'formulasMC' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.GenericFormulaAdder' ,
                  'declare'    : '',
                  'module'     : 'GenericFormulaAdder(\'data/formulasToAdd_MC.py\')' ,
                 },
   
  'formulasDATA' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.GenericFormulaAdder' ,
                  'declare'    : '',
                  'module'     : 'GenericFormulaAdder(\'data/formulasToAdd_DATA.py\')' ,
                 },

  'formulasFAKE' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'import'     : 'LatinoAnalysis.NanoGardener.modules.GenericFormulaAdder' ,
                  'declare'    : '',
                  'module'     : 'GenericFormulaAdder(\'data/formulasToAdd_FAKE.py\')' ,
                 },

# ------------------------------------ SKIMS : CUTS ONLY ----------------------------------------------------------

## ------- Fake Study:

  'fakeSel'    : {
                  'isChain'    : False ,
                  'do4MC'      : False  ,
                  'do4Data'    : True  ,
                  'selection'  : '"((MET_pt < 20 || PuppiMET_pt < 20) && mtw1 < 20)"' ,
                 },

  'fakeSelMC'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  , 
                  'selection'  : '"((MET_pt < 20 || PuppiMET_pt < 20) && mtw1 < 20)"' , 
                  'onlySample' : [
                                  #### DY
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  ####
                                  'WJetsToLNu','WJetsToLNu_HT100_200','WJetsToLNu_HT200_400','WJetsToLNu_HT400_600','WJetsToLNu_HT600_800',
                                  'WJetsToLNu_HT800_1200','WJetsToLNu_HT1200_2500','WJetsToLNu_HT2500_inf',
                                  ####
                                  'QCD_Pt-15to20_EMEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched','QCD_Pt-50to80_EMEnriched_ext1',
                                  'QCD_Pt-20toInf_MuEnrichedPt15','QCD_Pt-30toInf_DoubleEMEnriched','QCD_Pt-15to20_MuEnrichedPt5',
                                  ####
                                  'QCD_Pt_15to20_bcToE','QCD_Pt_20to30_bcToE','QCD_Pt_30to80_bcToE','QCD_Pt_80to170_bcToE',
                                  'QCD_Pt_170to250_bcToE','QCD_Pt_250toInf_bcToE',
                                  ####
                                  'TT','TTJets','TTTo2L2Nu',
                                 ] ,               
                 },

## ------- 2-Leptons: Loose / tightOR

  'l2loose'   : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  , 
                  'selection'  : '"(nLepton>=2)"' , 
                 },

  'l2tightOR2017' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'selection'  : '" (nLepton>=2 && Lepton_pt[0]>18 && Lepton_pt[1]>8 ) \
                                    && (    Lepton_isTightElectron_mvaFall17Iso_WP90[0] > 0.5        \
                                         || Lepton_isTightElectron_mvaFall17Iso_WP90_SS[0] > 0.5     \
                                         || Lepton_isTightMuon_cut_Tight_HWWW[0] > 0.5             ) \
                                    && (    Lepton_isTightElectron_mvaFall17Iso_WP90[1] > 0.5        \
                                         || Lepton_isTightElectron_mvaFall17Iso_WP90_SS[1] > 0.5     \
                                         || Lepton_isTightMuon_cut_Tight_HWWW[1] > 0.5             ) \
                                  "' , 
                 },

## ------- Analysis Skims:

# ------------------------------------ SPECIAL STEPS: HADD & UEPS -------------------------------------------------

## ------- HADD 

  'hadd'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'SizeMax'    : 5e9 ,
                 #'bigSamples' : ['DYJetsToLL_M-50','DY2JetsToLL','ZZTo2L2Q','DYJetsToLL_M-50-LO',
                 #                'DYJetsToLL_M-50-LO-ext1',
                 #                'WZTo2L2Q','TTToSemiLepton','TTToSemiLeptonic','TTTo2L2Nu_ext1','TTJetsDiLep-LO-ext1','TTTo2L2Nu',
                 #                'DYJetsToEE_Pow',
                 #                'DY1JetsToLL',
                 #                #'TTJets',
                 #               ],
               },

## ------- UEPS 

  'UEPS'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'onlySample' : [
                                    'GluGluHToWWTo2L2Nu_M125_CUETDown' , 'VBFHToWWTo2L2Nu_M125_CUETDown' , 'WWTo2L2Nu_CUETDown' ,
                                    'GluGluHToWWTo2L2Nu_M125_CUETUp'   , 'VBFHToWWTo2L2Nu_M125_CUETUp'   , 'WWTo2L2Nu_CUETUp'   ,
                                    'GluGluHToWWTo2L2NuHerwigPS_M125'  , 'VBFHToWWTo2L2NuHerwigPS_M125'  , 'WWTo2L2NuHerwigPS'  ,
                                    'GluGluHToWWTo2L2Nu_M125_herwigpp' ,
                                 ] ,
                  'cpMap' : {
                              'UEdo' : {
                                          'GluGluHToWWTo2L2Nu_M125_CUETDown' : ['GluGluHToWWTo2L2Nu_M125' ,'GluGluHToWWTo2L2NuPowheg_M125'],
                                          'VBFHToWWTo2L2Nu_M125_CUETDown'    : ['VBFHToWWTo2L2Nu_M125','VBFHToWWTo2L2Nu_alternative_M125']    ,
                                          'WWTo2L2Nu_CUETDown'               : ['WWTo2L2Nu'] ,
                                       },
                              'UEup' : {
                                          'GluGluHToWWTo2L2Nu_M125_CUETUp'   : ['GluGluHToWWTo2L2Nu_M125' ,'GluGluHToWWTo2L2NuPowheg_M125'],
                                          'VBFHToWWTo2L2Nu_M125_CUETUp'      : ['VBFHToWWTo2L2Nu_M125','VBFHToWWTo2L2Nu_alternative_M125']    ,
                                          'WWTo2L2Nu_CUETUp'                 : ['WWTo2L2Nu'] ,
                                       },
                              'PS'   : {
                                          'GluGluHToWWTo2L2NuHerwigPS_M125'  : ['GluGluHToWWTo2L2Nu_M125' ,'GluGluHToWWTo2L2NuPowheg_M125'],
                                          'GluGluHToWWTo2L2Nu_M125_herwigpp' : ['GluGluHToWWTo2L2Nu_M125' ,'GluGluHToWWTo2L2NuPowheg_M125'],
                                          'VBFHToWWTo2L2NuHerwigPS_M125'     : ['VBFHToWWTo2L2Nu_M125','VBFHToWWTo2L2Nu_alternative_M125'] ,
                                          'WWTo2L2NuHerwigPS'                : ['WWTo2L2Nu'] ,
                                       },
                            },
               },

} 


