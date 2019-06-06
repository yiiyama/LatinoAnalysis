import os
import sys
import collections
import array
import ROOT

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

class FakeSkimMaker(Module):
    def __init__(self, cmssw, skimNames, isMC, latinos=True):
        ROOT.gSystem.Load('libLatinoAnalysisNanoGardener.so')
        ROOT.gROOT.LoadMacro(os.path.dirname(__file__) + '/FakeSkimMaker.cc+')
        
        self.event = ROOT.fakeskim.Event(isMC)
        self.event.useLatinos(latinos)
        self.skims = [getattr(ROOT.fakeskim, sname)(self.event) for sname in skimNames]
        self.isMC = isMC

        self.electronMinPt = 13.
        self.electronMaxEta = 2.5
        self.muonMinPt = 10.
        self.muonMaxEta = 2.4

        # put the following in some data file
        if cmssw in ['Full2017v2']:
            self.setupFull2017v2()
        else:
            raise NotImplementedError('FakeSkimMaker not set up for production ' + cmssw)

    def beginJob(self):
        pass
        
    def endJob(self):
        pass
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        if self.event.usingLatinos():
            self.out.branch('Lepton_isBaseline', 'O', lenVar='nLepton')

        for skim in self.skims:
            self.out.branch('pass%s' % skim.getName(), 'O')
            skim.bookBranches(self.out._tree)

        self.event.bookBranches(self.out._tree)

        self.initReaders(inputTree) # initReaders must be called in beginFile

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def initReaders(self, tree):
        self.event.setTreeReader(tree._ttreereader)
        self._ttreereaderversion = tree._ttreereaderversion

    def analyze(self, event):
        skims = list(self.skims)

        if self.isMC:
            try:
                lheParts = Collection(event, 'LHEPart', lenVar='nLHEPart')
            except RuntimeError:
                # this sample does not have LHE
                pass
            else:
                nPart = collections.defaultdict(int)
                for part in lheParts:
                    nPart[part.pdgId] += 1
            
                iskim = 0
                while iskim != len(skims):
                    if skims[iskim].isGoodLHEEvent(nPart[11], nPart[-11], nPart[13], nPart[-13]):
                        iskim += 1
                    else:
                        skims.pop(iskim)
    
                if len(skims) == 0:
                    return False
        else:
            iskim = 0
            while iskim != len(skims):
                for trigger in skims[iskim].triggers():
                    try:
                        if event[str(trigger)]:
                            break
                    except RuntimeError:
                        pass
                else:
                    # no ORed trigger fired for this skim
                    skims.pop(iskim)
                    continue
    
                iskim += 1

            if len(skims) == 0:
                return False

        if self.event.usingLatinos():
            leptons = Collection(event, 'Lepton', lenVar='nLepton')
            electrons = Collection(event, 'Electron', lenVar='nElectron')
            muons = Collection(event, 'Muon', lenVar='nMuon')

            isBaseline = ROOT.vector('bool')(len(leptons), False)
            hasBaseline = False
    
            for il, lepton in enumerate(leptons):
                if lepton.electronIdx >= 0:
                    if self.isBaselineElectron(electrons[lepton.electronIdx], lepton):
                        isBaseline[il] = True
                        hasBaseline = True
                else:
                    if self.isBaselineMuon(muons[lepton.muonIdx], lepton):
                        isBaseline[il] = True
                        hasBaseline = True
    
            if not hasBaseline:
                return False
        else:
            isBaseline = None

        if event._tree._ttreereaderversion > self._ttreereaderversion:
            current = event._tree._ttreereader.GetCurrentEntry()
            event._tree._ttreereader.Restart()
            self.initReaders(event._tree)
            event._tree._ttreereader.SetEntry(current)

        if not self.event.selectLeptons(isBaseline):
            return False

        iskim = 0
        while iskim != len(skims):
            skim = skims[iskim]
            if skim.passSkim():
                iskim += 1
            else:
                skims.pop(iskim)

        if len(skims):
            return False

        self.event.setFlags()

        for skim in skims:
            skim.setWeights()

        if self.event.usingLatinos():
            self.out.fillBranch('Lepton_isBaseline', list(bool(x) for x in isBaseline))
        for skim in self.skims:
            self.out.fillBranch('pass%s' % skim.getName(), (skim in skims))

        return True

    def isBaselineElectron(self, electron, lepton):
        if electron.pt < self.electronMinPt or abs(electron.eta) > self.electronMaxEta:
            return False

        bitmap = electron.vidNestedWPBitmap

        for ith, thres in enumerate(self._eleBaselineRequirements):
            if thres == 0:
                continue
            if ((bitmap >> (ith * 3)) & 0x7) < thres:
                return False

        return True

    def isBaselineMuon(self, muon, lepton):
        return muon.pt > self.muonMinPt and abs(muon.eta) < self.muonMaxEta and lepton.isLoose

    def setupFull2017v2(self):
        self._eleBaselineRequirements = [
            0, # MinPtCut
            0, # GsfEleSCEtaMultiRangeCut
            1, # GsfEleDEtaInSeedCut
            1, # GsfEleDPhiInCut
            1, # GsfEleFull5x5SigmaIEtaIEtaCut
            1, # GsfEleHadronicOverEMEnergyScaledCut
            1, # GsfEleEInverseMinusPInverseCut
            1, # GsfEleRelPFIsoScaledCut
            0, # GsfEleConversionVetoCut
            0 # GsfEleMissingHitsCut
        ]

        self.event.leptons.addTightElectronWP('mvaFall17Iso_WP90')
        self.event.leptons.addTightElectronWP('mvaFall17Iso_WP90_SS')
        self.event.leptons.addTightMuonWP('cut_Tight_HWWW')

        cmssw_base = os.getenv('CMSSW_BASE')
        datadir = cmssw_base + '/src/LatinoAnalysis/NanoGardener/python/data'

        source = ROOT.TFile.Open(datadir + '/fake_decomposition/fall17_bhadron_flateff_60.root')
        btagMap = source.Get('varmap')
        self.event.setElectronBTagMap(btagMap)
        source.Close()

        if self.isMC:
            etabinning = [-2.5, -2.1, -1.6, -1.4, -0.8, 0., 0.8, 1.4, 1.6, 2.1, 2.5]
            ptbinning = [0., 10., 20., 30.] + [32. + x for x in range(7)] + [40., 45., 50., 60., 100., 200.]
            ele35WPTightGsfMap = ROOT.TH2D('ele35', '', len(ptbinning) - 1, array.array('d', ptbinning), len(etabinning) - 1, array.array('d', etabinning))
            ptbinning = [0., 10.] + [20. + x for x in range(7)] + [30. + 5. * x for x in range(5)] + [60., 100., 200.]
            ele23CaloIdLTrackIdLIsoVLMap = ROOT.TH2D('ele23', '', len(ptbinning) - 1, array.array('d', ptbinning), len(etabinning) - 1, array.array('d', etabinning))
            ptbinning = [0.] + [10. + x for x in range(6)] + [20. + 5. * x for x in range(7)] + [60., 100., 200.]
            ele12CaloIdLTrackIdLIsoVLMap = ROOT.TH2D('ele12', '', len(ptbinning) - 1, array.array('d', ptbinning), len(etabinning) - 1, array.array('d', etabinning))
    
            lumitotal = 0.
            for era, lumi in [('B', 4.823), ('CDE', 9.664 + 4.252 + 9.278), ('F', 13.54)]:
                lumitotal += lumi
                
                tmpmap = ele35WPTightGsfMap.Clone('tmpmap')
                with open(datadir + '/trigger/Full2017/Ele35_pt_eta_efficiency_withSys_Run2017' + era + '.txt') as source:
                    for line in source:
                        etamin, etamax, ptmin, ptmax, eff, _, _ = map(float, line.split())
                        ix = tmpmap.GetXaxis().FindFixBin((ptmin + ptmax) * 0.5)
                        iy = tmpmap.GetYaxis().FindFixBin((etamin + etamax) * 0.5)
                        tmpmap.SetBinContent(ix, iy, eff * lumi)
    
                ele35WPTightGsfMap.Add(tmpmap)
                tmpmap.Delete()
    
                tmpmap = ele23CaloIdLTrackIdLIsoVLMap.Clone('tmpmap')
                with open(datadir + '/trigger/Full2017/Ele23_Ele12_leg1_pt_eta_efficiency_withSys_Run2017' + era + '.txt') as source:
                    for line in source:
                        etamin, etamax, ptmin, ptmax, eff, _, _ = map(float, line.split())
                        ix = tmpmap.GetXaxis().FindFixBin((ptmin + ptmax) * 0.5)
                        iy = tmpmap.GetYaxis().FindFixBin((etamin + etamax) * 0.5)
                        tmpmap.SetBinContent(ix, iy, eff * lumi)
    
                ele23CaloIdLTrackIdLIsoVLMap.Add(tmpmap)
                tmpmap.Delete()
    
                tmpmap = ele12CaloIdLTrackIdLIsoVLMap.Clone('tmpmap')
                with open(datadir + '/trigger/Full2017/Ele23_Ele12_leg2_pt_eta_efficiency_withSys_Run2017' + era + '.txt') as source:
                    for line in source:
                        etamin, etamax, ptmin, ptmax, eff, _, _ = map(float, line.split())
                        ix = tmpmap.GetXaxis().FindFixBin((ptmin + ptmax) * 0.5)
                        iy = tmpmap.GetYaxis().FindFixBin((etamin + etamax) * 0.5)
                        tmpmap.SetBinContent(ix, iy, eff * lumi)
    
                ele12CaloIdLTrackIdLIsoVLMap.Add(tmpmap)
                tmpmap.Delete()
    
            ele35WPTightGsfMap.Scale(1. / lumitotal)
            self.event.setScaleFactorMap(ROOT.fakeskim.Event.Ele35_WPTight_Gsf, ele35WPTightGsfMap)
            ele23CaloIdLTrackIdLIsoVLMap.Scale(1. / lumitotal)
            self.event.setScaleFactorMap(ROOT.fakeskim.Event.Ele23_CaloIdL_TrackIdL_IsoVL, ele23CaloIdLTrackIdLIsoVLMap)
            ele12CaloIdLTrackIdLIsoVLMap.Scale(1. / lumitotal)
            self.event.setScaleFactorMap(ROOT.fakeskim.Event.Ele12_CaloIdL_TrackIdL_IsoVL, ele12CaloIdLTrackIdLIsoVLMap)
    
            etabinning = [-2.4, -2.1, -1.6, -1.2, -0.8, -0.3, -0.2, 0.2, 0.3, 0.8, 1.2, 1.6, 2.1, 2.4]
            ptbinning = [0., 10.] + [20. + x for x in range(7)] + [30. + 5. * x for x in range(5)] + [60., 100., 200.]
            mu23TrkIsoVVLMap = ROOT.TH2D('mu23', '', len(ptbinning) - 1, array.array('d', ptbinning), len(etabinning) - 1, array.array('d', etabinning))
            ptbinning = [0.] + [10. + x for x in range(6)] + [20. + 5. * x for x in range(7)] + [60., 100., 200.]
            mu12TrkIsoVVLMap = ROOT.TH2D('mu12', '', len(ptbinning) - 1, array.array('d', ptbinning), len(etabinning) - 1, array.array('d', etabinning))
    
            lumitotal = 0.
            for era, lumi in [('B', 4.823), ('CD', 9.664 + 4.252), ('E', 9.278), ('F', 13.54)]:
                lumitotal += lumi
                
                tmpmap = mu23TrkIsoVVLMap.Clone('tmpmap')
                with open(datadir + '/trigger/Full2017/Mu23_pt_eta_efficiency_withSys_Run2017' + era + '.txt') as source:
                    for line in source:
                        etamin, etamax, ptmin, ptmax, eff, _, _ = map(float, line.split())
                        ix = tmpmap.GetXaxis().FindFixBin((ptmin + ptmax) * 0.5)
                        iy = tmpmap.GetYaxis().FindFixBin((etamin + etamax) * 0.5)
                        tmpmap.SetBinContent(ix, iy, eff * lumi)
    
                mu23TrkIsoVVLMap.Add(tmpmap)
                tmpmap.Delete()
    
                tmpmap = mu12TrkIsoVVLMap.Clone('tmpmap')
                with open(datadir + '/trigger/Full2017/Mu12_pt_eta_efficiency_withSys_Run2017' + era + '.txt') as source:
                    for line in source:
                        etamin, etamax, ptmin, ptmax, eff, _, _ = map(float, line.split())
                        ix = tmpmap.GetXaxis().FindFixBin((ptmin + ptmax) * 0.5)
                        iy = tmpmap.GetYaxis().FindFixBin((etamin + etamax) * 0.5)
                        tmpmap.SetBinContent(ix, iy, eff * lumi)
    
                mu12TrkIsoVVLMap.Add(tmpmap)
                tmpmap.Delete()
    
            mu23TrkIsoVVLMap.Scale(1. / lumitotal)
            self.event.setScaleFactorMap(ROOT.fakeskim.Event.Mu23_TrkIsoVVL, mu23TrkIsoVVLMap)
            mu12TrkIsoVVLMap.Scale(1. / lumitotal)
            self.event.setScaleFactorMap(ROOT.fakeskim.Event.Mu12_TrkIsoVVL, mu12TrkIsoVVLMap)

            etabinning = [-2.5, -2., -1.566, -1.444, -0.8, 0., 0.8, 1.444, 1.566, 2., 2.5]
            ptbinning = [10., 20., 35., 50., 90., 150., 500.]
            eleTightIdMap = ROOT.TH2D('eleid', '', len(ptbinning) - 1, array.array('d', ptbinning), len(etabinning) - 1, array.array('d', etabinning))
    
            lumitotal = 0.
            for era, lumi in [('B', 4.823), ('C', 9.664), ('D', 4.252), ('E', 9.278), ('F', 13.54)]:
                lumitotal += lumi
                
                tmpmap = eleTightIdMap.Clone('tmpmap')
                with open(datadir + '/scale_factor/Full2017/egammaEffi_passingMVA94Xwp90isoHWW_run' + era + '.txt') as source:
                    for line in source:
                        entries = map(float, line.split())
                        etamin, etamax, ptmin, ptmax = entries[:4]
                        dataeff = entries[4]
                        mceff = entries[6]
                        ix = tmpmap.GetXaxis().FindFixBin((ptmin + ptmax) * 0.5)
                        iy = tmpmap.GetYaxis().FindFixBin((etamin + etamax) * 0.5)
                        tmpmap.SetBinContent(ix, iy, dataeff / mceff * lumi)
    
                eleTightIdMap.Add(tmpmap)
                tmpmap.Delete()
    
            eleTightIdMap.Scale(1. / lumitotal)
            self.event.setScaleFactorMap(ROOT.fakeskim.Event.ElectronTightId, eleTightIdMap)
