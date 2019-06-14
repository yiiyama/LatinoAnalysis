#include "LatinoAnalysis/NanoGardener/interface/DataFormats.h"

#include "TH2.h"
#include "TROOT.h"

#include <vector>
#include <utility>
#include <array>

namespace fakeskim {

  class PtEtaMap {
  public:
    PtEtaMap() {}
    ~PtEtaMap() {}
    void set(TH2 const*);
    double get(double pt, double eta) const;

  private:
    std::unique_ptr<TH2> map_{};
  };

  class Event {
  public:
    Event() {}
    Event(bool isMC) : isMC_(isMC) {}
    ~Event() {}

    //! Don't use precomputed Lepton branches
    void useLatinos(bool s) { usingLatinos_ = s; }
    void setTreeReader(TTreeReader*);

    void bookBranches(TTree*);
    void setFlags();

    bool isMC() const { return isMC_; }
    bool usingLatinos() const { return usingLatinos_; }

    nanoaod::LeptonCollection leptons{};
    nanoaod::ElectronCollection electrons{};
    nanoaod::MuonCollection muons{};
    nanoaod::JetCollection jets{};
    nanoaod::PhotonCollection photons{};
    nanoaod::MET met{};
    nanoaod::TrigObjCollection trigObjs{};
    nanoaod::GenPartCollection genParts{};

    UIntValueReaderPtr run{};
    UIntValueReaderPtr luminosityBlock{};
    ULong64ValueReaderPtr event{};

    enum ScaleFactor {
      ElectronTightId,
      Ele35_WPTight_Gsf,
      Ele23_CaloIdL_TrackIdL_IsoVL,
      Ele12_CaloIdL_TrackIdL_IsoVL,
      Mu23_TrkIsoVVL,
      Mu12_TrkIsoVVL,
      nScaleFactors
    };
    void setScaleFactorMap(ScaleFactor f, TH2 const* map) { scaleFactors_[f].set(map); }
    PtEtaMap const& getScaleFactorMap(ScaleFactor f) const { return scaleFactors_[f]; }

    void setElectronBTagMap(TH2 const* map) { electronBTagMap_.set(map); }

    bool selectLeptons(std::vector<bool> const* = nullptr);

    static constexpr unsigned NMAX{128};

    bool Electron_prompt[NMAX]{};
    bool Electron_conversion[NMAX]{};
    bool Electron_tau[NMAX]{};
    int Electron_hadronFlavour[NMAX]{};
    bool Electron_isBaseline[NMAX]{};
    bool Electron_isTight[NMAX]{};
    bool Electron_btag[NMAX]{};
    float Electron_btagDeepB[NMAX]{};
    float Electron_btagDeepC[NMAX]{};
    float Electron_jetPt[NMAX]{};
    float Electron_jetEta[NMAX]{};
    bool Muon_isBaseline[NMAX]{};
    bool Muon_isTight[NMAX]{};
    
  private:
    bool isMC_{};
    bool usingLatinos_{true};

    std::array<PtEtaMap, nScaleFactors> scaleFactors_{};
    PtEtaMap electronBTagMap_{};
  };

  class FakeSkim {
  public:
    FakeSkim(Event const& event) : event_(event) {}
    FakeSkim(FakeSkim const& orig) : event_(orig.event_), triggers_(orig.triggers_) {}
    virtual ~FakeSkim() {}
    virtual char const* getName() const = 0;
    virtual void bookBranches(TTree*) = 0;
    std::vector<TString> const& triggers() const { return triggers_; }
    virtual bool isGoodLHEEvent(unsigned nElectron, unsigned nPositron, unsigned nMuon, unsigned nAntimuon) const = 0;
    virtual bool passSkim() = 0;
    virtual void setWeights() {}
  protected:
    Event const& event_;
    std::vector<TString> triggers_;
  };

  class OppositeSignDielectronSkim : public FakeSkim {
  public:
    OppositeSignDielectronSkim(Event const& event);
    OppositeSignDielectronSkim(OppositeSignDielectronSkim const&);
    ~OppositeSignDielectronSkim() {}
    char const* getName() const override { return "OppositeSignDielectronSkim"; }
    void bookBranches(TTree* outTree) override;
    bool isGoodLHEEvent(unsigned nElectron, unsigned nPositron, unsigned nMuon, unsigned nAntimuon) const override {
      return nElectron == 1 && nPositron == 1 && nMuon == 0 && nAntimuon == 0;
    }
    bool passSkim() override;

  private:
    bool Electron_isTag[Event::NMAX]{};
    bool Electron_isProbe[Event::NMAX]{};
    bool Electron_isCaloIdLTrackIdLIsoVL[Event::NMAX]{};
    float Electron_mee[Event::NMAX]{};
    float Electron_tagWeight[Event::NMAX]{};
    unsigned Electron_njet[Event::NMAX]{};
    unsigned Electron_nbjet[Event::NMAX]{};
  };

  class SameSignDielectronSkim : public FakeSkim {
  public:
    SameSignDielectronSkim(Event const& event);
    SameSignDielectronSkim(SameSignDielectronSkim const&);
    ~SameSignDielectronSkim() {}
    char const* getName() const override { return "SameSignDielectronSkim"; }
    void bookBranches(TTree* outTree) override;
    bool isGoodLHEEvent(unsigned nElectron, unsigned nPositron, unsigned nMuon, unsigned nAntimuon) const override {
      return nElectron + nPositron >= 1;
    }
    bool passSkim() override;
    void setWeights() override;

  private:
    std::array<unsigned, 2> indices{};
    unsigned njet{};
    unsigned nbjet{};
    float skimWeight{};
  };

  class DimuonElectronSkim : public FakeSkim {
  public:
    DimuonElectronSkim(Event const& event);
    DimuonElectronSkim(DimuonElectronSkim const&);
    ~DimuonElectronSkim() {}
    char const* getName() const override { return "DimuonElectronSkim"; }
    void bookBranches(TTree* outTree) override;
    bool isGoodLHEEvent(unsigned nElectron, unsigned nPositron, unsigned nMuon, unsigned nAntimuon) const override {
      return (nMuon == 1 || nAntimuon == 1);
    }
    bool passSkim() override;

  private:
    std::array<unsigned, 2> tightMu{};
    float mmm{};
    float Electron_minDRMu[Event::NMAX]{};
    float Electron_isCaloIdLTrackIdLIsoVL[Event::NMAX]{};
    float Electron_mmme[Event::NMAX]{};
    unsigned Electron_njet[Event::NMAX]{};
    unsigned Electron_nbjet[Event::NMAX]{};
    float skimWeight{};
  };

  class DimuonPhotonSkim : public FakeSkim {
  public:
    DimuonPhotonSkim(Event const& event);
    DimuonPhotonSkim(DimuonPhotonSkim const&);
    ~DimuonPhotonSkim() {}
    char const* getName() const override { return "DimuonPhotonSkim"; }
    void bookBranches(TTree* outTree) override;
    bool isGoodLHEEvent(unsigned nElectron, unsigned nPositron, unsigned nMuon, unsigned nAntimuon) const override {
      return (nMuon == 1 || nAntimuon == 1);
    }
    bool passSkim() override;

  private:
    std::array<unsigned, 2> tightMu{};
    float mmm{};
    float Photon_minDRMu[Event::NMAX]{};
    float Photon_mmmg[Event::NMAX]{};
    unsigned Photon_njet[Event::NMAX]{};
    unsigned Photon_nbjet[Event::NMAX]{};
    float skimWeight{};
  };

  class SameSignMuonElectronSkim : public FakeSkim {
  public:
    SameSignMuonElectronSkim(Event const& event);
    SameSignMuonElectronSkim(SameSignMuonElectronSkim const&);
    ~SameSignMuonElectronSkim() {}
    char const* getName() const override { return "SameSignMuonElectronSkim"; }
    void bookBranches(TTree* outTree) override;
    bool isGoodLHEEvent(unsigned nElectron, unsigned nPositron, unsigned nMuon, unsigned nAntimuon) const override {
      return nElectron + nPositron + nMuon + nAntimuon >= 1;
    }
    bool passSkim() override;

  private:
    unsigned tightMu{};
    unsigned njet{};
    unsigned nbjet{};
    float skimWeight{};
  };

  class JetElectronSkim : public FakeSkim {
  public:
    JetElectronSkim(Event const& event);
    JetElectronSkim(JetElectronSkim const&);
    ~JetElectronSkim() {}
    char const* getName() const override { return "JetElectronSkim"; }
    void bookBranches(TTree* outTree) override;
    bool isGoodLHEEvent(unsigned nElectron, unsigned nPositron, unsigned nMuon, unsigned nAntimuon) const override {
      return true;
    }
    bool passSkim() override;

  private:
    float Electron_minDPhiJet[Event::NMAX]{};
    float Electron_ptOppJet[Event::NMAX]{};
    unsigned Electron_njet[Event::NMAX]{};
    unsigned Electron_nbjet[Event::NMAX]{};
  };

  class ElectronSkim : public FakeSkim {
  public:
    ElectronSkim(Event const& event);
    ElectronSkim(ElectronSkim const&);
    ~ElectronSkim() {}
    char const* getName() const override { return "ElectronSkim"; }
    void bookBranches(TTree* outTree) override;
    bool isGoodLHEEvent(unsigned nElectron, unsigned nPositron, unsigned nMuon, unsigned nAntimuon) const override {
      return true;
    }
    bool passSkim() override;

  private:
    float Electron_isCaloIdLTrackIdLIsoVL[Event::NMAX]{};
  };

}

void
fakeskim::PtEtaMap::set(TH2 const* map)
{
  gROOT->cd();
  map_.reset(static_cast<TH2*>(map->Clone()));
  map_->SetDirectory(nullptr);
}

double
fakeskim::PtEtaMap::get(double pt, double eta) const
{
  int ix(map_->GetXaxis()->FindFixBin(pt));
  if (ix == 0)
    ix = 1;
  else if (ix >= map_->GetNbinsX())
    ix = map_->GetNbinsX();
  int iy(map_->GetYaxis()->FindFixBin(eta));
  if (iy == 0)
    iy = 1;
  else if (iy >= map_->GetNbinsY())
    iy = map_->GetNbinsY();

  return map_->GetBinContent(ix, iy);
}

void
fakeskim::Event::setTreeReader(TTreeReader* reader)
{
  //hlt.setTreeReader(_reader, "HLT");
  if (usingLatinos_)
    leptons.setTreeReader(reader, "Lepton");
  electrons.setTreeReader(reader, "Electron");
  muons.setTreeReader(reader, "Muon");
  jets.setTreeReader(reader, "Jet");
  photons.setTreeReader(reader, "Photon");
  trigObjs.setTreeReader(reader, "TrigObj");
  met.setTreeReader(reader, "MET");
  if (isMC_) {
    genParts.setTreeReader(reader, "GenPart");
    //lheParts.setTreeReader(reader, "LHEPart");
  }

  run = std::make_unique<UIntValueReader>(*reader, "run");
  luminosityBlock = std::make_unique<UIntValueReader>(*reader, "luminosityBlock");
  event = std::make_unique<ULong64ValueReader>(*reader, "event");
  // if (isMC_)
  //   genWeight = std::make_unique<FloatValueReader>(*reader, "genWeight");
}

void
fakeskim::Event::bookBranches(TTree* outTree)
{
  outTree->Branch("Electron_prompt", Electron_prompt, "Electron_prompt[nElectron]/O");
  outTree->Branch("Electron_conversion", Electron_conversion, "Electron_conversion[nElectron]/O");
  outTree->Branch("Electron_tau", Electron_tau, "Electron_tau[nElectron]/O");
  outTree->Branch("Electron_hadronFlavour", Electron_hadronFlavour, "Electron_hadronFlavour[nElectron]/I");
  outTree->Branch("Electron_isBaseline", Electron_isBaseline, "Electron_isBaseline[nElectron]/O");
  outTree->Branch("Electron_isTight", Electron_isTight, "Electron_isTight[nElectron]/O");
  outTree->Branch("Electron_btag", Electron_btag, "Electron_btag[nElectron]/O");
  outTree->Branch("Electron_btagDeepB", Electron_btagDeepB, "Electron_btagDeepB[nElectron]/F");
  outTree->Branch("Electron_btagDeepC", Electron_btagDeepC, "Electron_btagDeepC[nElectron]/F");
  outTree->Branch("Electron_jetPt", Electron_jetPt, "Electron_jetPt[nElectron]/F");
  outTree->Branch("Electron_jetEta", Electron_jetEta, "Electron_jetEta[nElectron]/F");
  outTree->Branch("Muon_isBaseline", Muon_isBaseline, "Muon_isBaseline[nMuon]/O");
  outTree->Branch("Muon_isTight", Muon_isTight, "Muon_isTight[nMuon]/O");
}

void
fakeskim::Event::setFlags()
{
  unsigned nE(electrons.size());

  if (isMC_) {
    std::fill_n(Electron_prompt, nE, false);
    std::fill_n(Electron_conversion, nE, false);
    std::fill_n(Electron_tau, nE, false);
    std::fill_n(Electron_hadronFlavour, nE, 0);

    std::vector<unsigned> ielectrons;
    std::vector<unsigned> iphotons;
    std::vector<unsigned> itaus;
    std::vector<unsigned> ibhadrons;
    std::vector<unsigned> ichadrons;

    unsigned nG(genParts.size());

    for (unsigned iG(0); iG != nG; ++iG) {
      unsigned absId(std::abs(genParts.pdgId(iG)));

      if (((absId / 100) % 10) == 5 || ((absId / 1000) % 10) == 5)
        ibhadrons.push_back(iG);
      else if (((absId / 100) % 10) == 4 || ((absId / 1000) % 10) == 4)
        ichadrons.push_back(iG);
      else if (absId == 15)
        itaus.push_back(iG);

      if (genParts.status(iG) != 1 || !genParts.checkStatus(iG, nanoaod::GenPartCollection::kIsPrompt, nanoaod::GenPartCollection::kIsPromptTauDecayProduct))
        continue;

      if (absId == 11)
        ielectrons.push_back(iG);
      else if (absId == 22 && genParts.checkStatus(iG, nanoaod::GenPartCollection::kIsPrompt) && genParts.pt(iG) > 10.)
        iphotons.push_back(iG);
    }

    for (unsigned iE(0); iE != nE; ++iE) {
      if (!Electron_isBaseline[iE])
        continue;

      Electron_prompt[iE] = electrons.deltaR2Match(iE, genParts, ielectrons, 0.01);
      Electron_conversion[iE] = electrons.deltaR2PtMatch(iE, genParts, iphotons, 0.0225, 0.5);
      Electron_tau[iE] = electrons.deltaR2Match(iE, genParts, itaus, 0.01);
      int jetIdx(electrons.jetIdx(iE));
      if (jetIdx >= 0)
        Electron_hadronFlavour[iE] = jets.hadronFlavour(jetIdx);
      else {
        // only for electrons with no matchin jet
        if (electrons.deltaR2Match(iE, genParts, ibhadrons, 0.16))
          Electron_hadronFlavour[iE] = 5;
        else if (electrons.deltaR2Match(iE, genParts, ichadrons, 0.16))
          Electron_hadronFlavour[iE] = 4;
        else
          Electron_hadronFlavour[iE] = 0;
      }
    }
  }
  
  std::fill_n(Electron_btag, nE, false);
  std::fill_n(Electron_btagDeepB, nE, -2.);
  std::fill_n(Electron_btagDeepC, nE, -1.);
  std::fill_n(Electron_jetPt, nE, 0.);
  std::fill_n(Electron_jetEta, nE, 0.);

  for (unsigned iE(0); iE != nE; ++iE) {
    if (!Electron_isBaseline[iE])
      continue;

    int jetIdx(electrons.jetIdx(iE));
    if (jetIdx < 0)
      continue;

    Electron_btagDeepB[iE] = jets.btagDeepB(jetIdx);
    Electron_btagDeepC[iE] = jets.btagDeepC(jetIdx);
    Electron_btag[iE] = Electron_btagDeepB[iE] > electronBTagMap_.get(electrons.pt(iE), std::abs(electrons.eta(iE)));
    Electron_jetPt[iE] = jets.pt(jetIdx);
    Electron_jetEta[iE] = jets.eta(jetIdx);
  }
}

bool
fakeskim::Event::selectLeptons(std::vector<bool> const* flags/* = nullptr*/)
{
  unsigned nE(electrons.size());
  unsigned nM(muons.size());
  std::fill_n(Electron_isBaseline, nE, false);
  std::fill_n(Muon_isBaseline, nM, false);
  std::fill_n(Electron_isTight, nE, false);
  std::fill_n(Muon_isTight, nM, false);

  bool hasBaseline(false);

  if (usingLatinos_) {
    unsigned nL(leptons.size());
    for (unsigned iL(0); iL != nL; ++iL) {
      if (!(*flags)[iL])
        continue;

      hasBaseline = true;

      if (leptons.electronIdx(iL) >= 0) {
        Electron_isBaseline[leptons.electronIdx(iL)] = true;
        Electron_isTight[leptons.electronIdx(iL)] = (leptons.isLoose(iL) != 0 && leptons.isTight(iL));
      }
      else {
        Muon_isBaseline[leptons.muonIdx(iL)] = true;
        Muon_isTight[leptons.muonIdx(iL)] = (leptons.isLoose(iL) != 0 && leptons.isTight(iL));
      }
    }
  }
  else {
    unsigned interestingCuts(
      (1 << nanoaod::ElectronCollection::MinPtCut) +
      (1 << nanoaod::ElectronCollection::GsfEleSCEtaMultiRangeCut) +
      (1 << nanoaod::ElectronCollection::GsfEleDEtaInSeedCut) +
      (1 << nanoaod::ElectronCollection::GsfEleDPhiInCut) +
      (1 << nanoaod::ElectronCollection::GsfEleFull5x5SigmaIEtaIEtaCut) +
      (1 << nanoaod::ElectronCollection::GsfEleHadronicOverEMEnergyScaledCut) +
      (1 << nanoaod::ElectronCollection::GsfEleEInverseMinusPInverseCut) +
      (1 << nanoaod::ElectronCollection::GsfEleRelPFIsoScaledCut)// +
      //(1 << nanoaod::ElectronCollection::GsfEleConversionVetoCut) +
      //(1 << nanoaod::ElectronCollection::GsfEleMissingHitsCut)
    );

    for (unsigned iE(0); iE != nE; ++iE) {
      double absEta(std::abs(electrons.eta(iE)));

      if (electrons.pt(iE) < 13. || absEta > 2.5)
        continue;

      if (!electrons.passCuts(iE, nanoaod::ElectronCollection::Veto, interestingCuts))
        continue;

      hasBaseline = true;

      Electron_isBaseline[iE] = true;

      if (electrons.cutBased(iE) == 0 || !electrons.mvaFall17Iso_WP90(iE) || !electrons.convVeto(iE))
        continue;

      if (absEta < 1.479) {
        if (std::abs(electrons.dxy(iE)) > 0.05 || std::abs(electrons.dz(iE)) > 0.1)
          continue;
      }
      else if (absEta < 2.5) {
        if (electrons.sieie(iE) > 0.03 || electrons.eInvMinusPInv(iE) > 0.014)
          continue;
        if (std::abs(electrons.dxy(iE)) > 0.1 || std::abs(electrons.dz(iE)) > 0.2)
          continue;
      }

      Electron_isTight[iE] = true;
    }

    for (unsigned iM(0); iM != nM; ++iM) {
      if (!muons.tightId(iM) || std::abs(muons.dz(iM)) > 0.1 || muons.pfRelIso04_all(iM) > 0.4)
        continue;

      if (std::abs(muons.eta(iM)) > 2.4)
        continue;

      double pt(muons.pt(iM));
      if (pt < 10.)
        continue;
      else if (pt < 20.) {
        if (std::abs(muons.dxy(iM)) > 0.01)
          continue;
      }
      else {
        if (std::abs(muons.dxy(iM)) > 0.02)
          continue;
      }

      hasBaseline = true;

      Muon_isBaseline[iM] = true;

      if (muons.pfRelIso04_all(iM) > 0.15)
        continue;

      Muon_isTight[iM] = true;
    }
  }

  return hasBaseline;
}

fakeskim::OppositeSignDielectronSkim::OppositeSignDielectronSkim(Event const& event) :
  FakeSkim(event)
{
  triggers_.push_back("HLT_Ele35_WPTight_Gsf");
}

fakeskim::OppositeSignDielectronSkim::OppositeSignDielectronSkim(OppositeSignDielectronSkim const& orig) :
  FakeSkim(orig)
{
  triggers_.push_back("HLT_Ele35_WPTight_Gsf");
}

void
fakeskim::OppositeSignDielectronSkim::bookBranches(TTree* outTree)
{
  outTree->Branch("Electron_OS2E_isTag", Electron_isTag, "Electron_OS2E_isTag[nElectron]/O");
  outTree->Branch("Electron_OS2E_isProbe", Electron_isProbe, "Electron_OS2E_isProbe[nElectron]/O");
  outTree->Branch("Electron_OS2E_isCaloIdLTrackIdLIsoVL", Electron_isCaloIdLTrackIdLIsoVL, "Electron_OS2E_isCaloIdLTrackIdLIsoVL[nElectron]/O");
  outTree->Branch("Electron_OS2E_mee", Electron_mee, "Electron_OS2E_mee[nElectron]/F");
  outTree->Branch("Electron_OS2E_tagWeight", Electron_tagWeight, "Electron_OS2E_tagWeight[nElectron]/F");
  outTree->Branch("Electron_OS2E_njet", Electron_njet, "Electron_OS2E_njet[nElectron]/i");
  outTree->Branch("Electron_OS2E_nbjet", Electron_nbjet, "Electron_OS2E_nbjet[nElectron]/i");
}

bool
fakeskim::OppositeSignDielectronSkim::passSkim()
{
  auto& electrons(event_.electrons);
  unsigned nE(electrons.size());

  unsigned nBaselineElectrons(std::count(event_.Electron_isBaseline, event_.Electron_isBaseline + nE, true));
  unsigned nTightElectrons(std::count(event_.Electron_isTight, event_.Electron_isTight + nE, true));

  if (nBaselineElectrons != 2 || nTightElectrons == 0)
    return false;

  auto& trigObjs(event_.trigObjs);
  std::vector<unsigned> ele35Objs;
  std::vector<unsigned> ele23Objs;
  std::vector<unsigned> ele12Objs;
  unsigned nO(trigObjs.size());
  for (unsigned iO(0); iO != nO; ++iO) {
    if (trigObjs.id(iO) != 11)
      continue;
    int filterBits(trigObjs.filterBits(iO));
    double pt(trigObjs.pt(iO));
    if ((filterBits & 2) != 0 && pt > 35.)
      ele35Objs.push_back(iO);
    if ((filterBits & 1) != 0) {
      if (pt > 23.)
        ele23Objs.push_back(iO);
      if (pt > 12.)
        ele12Objs.push_back(iO);
    }
  }

  auto& jets(event_.jets);
  unsigned nJ(jets.size());
  unsigned njet(0);
  unsigned nbjet(0);
  for (unsigned iJ(0); iJ != nJ; ++iJ) {
    if (jets.pt(iJ) < 30. || std::abs(jets.eta(iJ)) > 4.7)
      continue;

    ++njet;

    if (jets.btagDeepB(iJ) > 0.4941)
      ++nbjet;
  }

  bool hasTagAndProbePair(false);

  std::fill_n(Electron_isTag, nE, false);
  std::fill_n(Electron_isProbe, nE, false);
  std::fill_n(Electron_isCaloIdLTrackIdLIsoVL, nE, false);
  std::fill_n(Electron_mee, nE, 0.);
  std::fill_n(Electron_tagWeight, nE, 0.);
  std::fill_n(Electron_njet, nE, 0);
  std::fill_n(Electron_nbjet, nE, 0);
  for (unsigned iTag(0); iTag != nE; ++iTag) {
    if (!event_.Electron_isTight[iTag])
      continue;

    double tagPt(electrons.pt(iTag));

    if (event_.isMC())
      Electron_isTag[iTag] = tagPt > 40.;
    else
      Electron_isTag[iTag] = tagPt > 40. && electrons.deltaR2Match(iTag, trigObjs, ele35Objs, 0.1);

    if (!Electron_isTag[iTag])
      continue;

    unsigned tagJetIdx(electrons.jetIdx(iTag));
    bool tagIsJet(false);
    bool tagIsBjet(false);
    if (tagJetIdx < nJ && jets.pt(tagJetIdx) > 30. && std::abs(jets.eta(tagJetIdx)) < 4.7) {
      tagIsJet = true;
      if (jets.btagDeepB(tagJetIdx) > 0.4941)
        tagIsBjet = true;
    }

    int tagPdgId(electrons.pdgId(iTag));
    
    for (unsigned iProbe(0); iProbe != nE; ++iProbe) {
      if (iProbe == iTag)
        continue;

      if (Electron_isProbe[iProbe])
        continue;

      if (!event_.Electron_isBaseline[iProbe])
        continue;

      if (electrons.pdgId(iProbe) * tagPdgId != -121)
        continue;

      hasTagAndProbePair = true;

      double probePt(electrons.pt(iProbe));

      Electron_isProbe[iProbe] = true;

      // We normalize the Z peak by all baseline electrons, but are interested in trigger + tag combined efficiencies.
      // Therefore we flag the electrons matching the trigger but won't continue / return false upon mismatch
      if (event_.isMC())
        Electron_isCaloIdLTrackIdLIsoVL[iProbe] = true;
      else
        Electron_isCaloIdLTrackIdLIsoVL[iProbe] = \
          electrons.deltaR2Match(iTag, trigObjs, ele23Objs, 0.1) ||
          electrons.deltaR2Match(iTag, trigObjs, ele12Objs, 0.1);

      Electron_mee[iProbe] = (electrons.p4(iTag) + electrons.p4(iProbe)).M();
      if (event_.isMC()) {
        Electron_tagWeight[iProbe] = \
          event_.getScaleFactorMap(Event::Ele35_WPTight_Gsf).get(tagPt, electrons.eta(iTag)) * \
          event_.getScaleFactorMap(Event::ElectronTightId).get(tagPt, electrons.etaSC(iTag));
      }

      Electron_njet[iProbe] = njet - (tagIsJet ? 1 : 0);

      unsigned probeJetIdx(electrons.jetIdx(iProbe));
      if (probeJetIdx < nJ && probeJetIdx != tagJetIdx && jets.pt(probeJetIdx) > 30. && std::abs(jets.eta(probeJetIdx)) < 4.7) {
        --Electron_njet[iProbe];
        if (jets.btagDeepB(probeJetIdx) > 0.4941)
          --Electron_nbjet[iProbe];
      }
    }
  }

  return hasTagAndProbePair;
}

fakeskim::SameSignDielectronSkim::SameSignDielectronSkim(Event const& event) :
  FakeSkim(event)
{
  triggers_.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
}

fakeskim::SameSignDielectronSkim::SameSignDielectronSkim(SameSignDielectronSkim const& orig) :
  FakeSkim(orig)
{
  triggers_.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
}

void
fakeskim::SameSignDielectronSkim::bookBranches(TTree* outTree)
{
  outTree->Branch("leadE_SS2E", indices.data(), "leadE_SS2E/i");
  outTree->Branch("trailE_SS2E", indices.data() + 1, "trailE_SS2E/i");
  outTree->Branch("njet_SS2E", &njet, "njet_SS2E/i");
  outTree->Branch("nbjet_SS2E", &nbjet, "nbjet_SS2E/i");
  outTree->Branch("skimWeight_SS2E", &skimWeight, "skimWeight_SS2E/F");
}

bool
fakeskim::SameSignDielectronSkim::passSkim()
{
  // TrigObj electrons require *PixelMatchFilter in the filter name whereas the dielectron filters have
  // name patterns *PixelMatchLeg{1,2}Filter (at least in 9_4_9, have not investigated later releases)

  // auto& trigObjs(event_.trigObjs);
  // std::vector<unsigned> ele23Objs;
  // std::vector<unsigned> ele12Objs;
  // if (!event_.isMC()) {
  //   unsigned nO(trigObjs.size());
  //   for (unsigned iO(0); iO != nO; ++iO) {
  //     if (trigObjs.id(iO) != 11)
  //       continue;
  //     int filterBits(trigObjs.filterBits(iO));
  //     double pt(trigObjs.pt(iO));
  //     if ((filterBits & 1) != 0) {
  //       if (pt > 23.)
  //         ele23Objs.push_back(iO);
  //       if (pt > 12.)
  //         ele12Objs.push_back(iO);
  //     }
  //   }
  // }

  auto& electrons(event_.electrons);
  unsigned nE(electrons.size());

  indices = {nE, nE};
  int leadPdgId(0);

  for (unsigned iE(0); iE != nE; ++iE) {
    if (!event_.Electron_isBaseline[iE])
      continue;

    // We cannot require trigger object matching (see above)
    // // Here we strictly require that there be exactly two baseline electrons, which are same sign, and match the trigger objects
    // if (!event_.isMC()) {
    //   if (!electrons.deltaR2Match(iE, trigObjs, ele23Objs, 0.1) && !electrons.deltaR2Match(iE, trigObjs, ele12Objs, 0.1))
    //     return false;
    // }

    if (indices[0] == nE)
      indices[0] = iE;
    else {
      if (indices[1] != nE || electrons.charge(iE) != electrons.charge(indices[0]))
        return false;

      indices[1] = iE;
    }
  }

  if (indices[1] == nE)
    return false;

  auto& jets(event_.jets);
  unsigned nJ(jets.size());
  njet = 0;
  nbjet = 0;
  unsigned jetIdx0(electrons.jetIdx(indices[0]));
  unsigned jetIdx1(electrons.jetIdx(indices[1]));
  for (unsigned iJ(0); iJ != nJ; ++iJ) {
    if (iJ == jetIdx0 || iJ == jetIdx1)
      continue;

    if (jets.pt(iJ) < 30. || std::abs(jets.eta(iJ)) > 4.7)
      continue;

    ++njet;

    if (jets.btagDeepB(iJ) > 0.4941)
      ++nbjet;
  }

  return true;
}

void
fakeskim::SameSignDielectronSkim::setWeights()
{
  if (event_.isMC()) {
    auto& electrons(event_.electrons);
    unsigned nE(electrons.size());

    skimWeight = 1.;
    
    for (unsigned iE(0); iE != nE; ++iE) {
      if (!event_.Electron_prompt[iE] || !event_.Electron_isTight[iE])
        continue;

      double pt(electrons.pt(iE));
      double eta(electrons.eta(iE));

      double nonFiring(1. - event_.getScaleFactorMap(Event::Ele23_CaloIdL_TrackIdL_IsoVL).get(pt, eta));
      nonFiring *= 1. - event_.getScaleFactorMap(Event::Ele12_CaloIdL_TrackIdL_IsoVL).get(pt, eta);
      skimWeight *= 1. - nonFiring;
      
      skimWeight *= event_.getScaleFactorMap(Event::ElectronTightId).get(pt, electrons.etaSC(iE));
    }
  }
}

fakeskim::DimuonElectronSkim::DimuonElectronSkim(Event const& event) :
  FakeSkim(event)
{
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  triggers_.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
}

fakeskim::DimuonElectronSkim::DimuonElectronSkim(DimuonElectronSkim const& orig) :
  FakeSkim(orig)
{
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  triggers_.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
}

void
fakeskim::DimuonElectronSkim::bookBranches(TTree* outTree)
{
  outTree->Branch("tightMu0_2ME", tightMu.data(), "tightMu0_2ME/i");
  outTree->Branch("tightMu1_2ME", tightMu.data() + 1, "tightMu1_2ME/i");
  outTree->Branch("mmm_2ME", &mmm, "mmm_2ME/F");
  outTree->Branch("Electron_2ME_isCaloIdLTrackIdLIsoVL", Electron_isCaloIdLTrackIdLIsoVL, "Electron_2ME_isCaloIdLTrackIdLIsoVL[nElectron]/O");
  outTree->Branch("Electron_2ME_minDRMu", Electron_minDRMu, "Electron_2ME_minDRMu[nElectron]/F");
  outTree->Branch("Electron_2ME_mmme", Electron_mmme, "Electron_2ME_mmme[nElectron]/F");
  outTree->Branch("Electron_2ME_njet", Electron_njet, "Electron_2ME_njet[nElectron]/i");
  outTree->Branch("Electron_2ME_nbjet", Electron_nbjet, "Electron_2ME_nbjet[nElectron]/i");
  outTree->Branch("skimWeight_2ME", &skimWeight, "skimWeight_2ME/F");
}

bool
fakeskim::DimuonElectronSkim::passSkim()
{
  auto& electrons(event_.electrons);
  unsigned nE(electrons.size());

  if (std::count(event_.Electron_isBaseline, event_.Electron_isBaseline + nE, true) == 0)
    return false;

  auto& muons(event_.muons);

  unsigned nM(muons.size());
  
  tightMu = {nM, nM};

  for (unsigned iM(0); iM != nM; ++iM) {
    if (!event_.Muon_isTight[iM])
      continue;

    if (tightMu[0] == nM)
      tightMu[0] = iM;
    else {
      if (tightMu[1] != nM || muons.charge(iM) == muons.charge(tightMu[0]))
        return false;
      
      tightMu[1] = iM;
    }
  }

  if (tightMu[1] == nM)
    return false;

  auto pmm(muons.p4(tightMu[0]) + muons.p4(tightMu[1]));
  mmm = pmm.M();
  if (mmm < 20.)
    return false;

  auto& trigObjs(event_.trigObjs);
  std::vector<unsigned> ele23Objs;
  std::vector<unsigned> ele12Objs;
  if (!event_.isMC()) {
    unsigned nO(trigObjs.size());
    for (unsigned iO(0); iO != nO; ++iO) {
      if (trigObjs.id(iO) != 11)
        continue;
      int filterBits(trigObjs.filterBits(iO));
      double pt(trigObjs.pt(iO));
      if ((filterBits & 1) != 0) {
        if (pt > 23.)
          ele23Objs.push_back(iO);
        if (pt > 12.)
          ele12Objs.push_back(iO);
      }
    }
  }

  std::fill_n(Electron_minDRMu, nE, 0.);
  std::fill_n(Electron_isCaloIdLTrackIdLIsoVL, nE, false);
  std::fill_n(Electron_mmme, nE, 0.);
  std::fill_n(Electron_njet, nE, 0);
  std::fill_n(Electron_nbjet, nE, 0);

  bool hasElectron(false);

  for (unsigned iE(0); iE != nE; ++iE) {
    if (!event_.Electron_isBaseline[iE])
      continue;

    // Here we are interested in tag + trigger combined efficiencies but can allow multiple baseline electrons.
    // The normalization for conversion scale factor efficiency comes from mumugamma, so we don't need to tag & save
    // but can simply continue in case of mismatch.
    if (!event_.isMC()) {
      if (!electrons.deltaR2Match(iE, trigObjs, ele23Objs, 0.1) && !electrons.deltaR2Match(iE, trigObjs, ele12Objs, 0.1))
        continue;
    }

    Electron_isCaloIdLTrackIdLIsoVL[iE] = true;
    hasElectron = true;

    Electron_minDRMu[iE] = std::sqrt(std::min(electrons.deltaR2(iE, muons, tightMu[0]), electrons.deltaR2(iE, muons, tightMu[1])));
    Electron_mmme[iE] = (pmm + electrons.p4(iE)).M();
  }

  if (!hasElectron)
    return false;

  auto& jets(event_.jets);
  unsigned nJ(jets.size());
  unsigned njet(0);
  unsigned nbjet(0);
  unsigned jetIdxMu0(muons.jetIdx(tightMu[0]));
  unsigned jetIdxMu1(muons.jetIdx(tightMu[1]));
  for (unsigned iJ(0); iJ != nJ; ++iJ) {
    if (iJ == jetIdxMu0 || iJ == jetIdxMu1)
      continue;

    if (jets.pt(iJ) < 30. || std::abs(jets.eta(iJ)) > 4.7)
      continue;

    ++njet;

    if (jets.btagDeepB(iJ) > 0.4941)
      ++nbjet;
  }

  for (unsigned iE(0); iE != nE; ++iE) {
    if (!event_.Electron_isBaseline[iE])
      continue;

    Electron_njet[iE] = njet;
    Electron_nbjet[iE] = nbjet;

    unsigned jetIdx(electrons.jetIdx(iE));
    if (jetIdx < nJ && jetIdx != jetIdxMu0 && jetIdx != jetIdxMu1 && jets.pt(jetIdx) > 30. && std::abs(jets.eta(jetIdx)) < 4.7) {
      --Electron_njet[iE];
      if (jets.btagDeepB(jetIdx) > 0.4941)
        --Electron_nbjet[iE];
    }
  }

  if (event_.isMC()) {
    double pt0(muons.pt(tightMu[0]));
    double eta0(muons.eta(tightMu[0]));
    double pt1(muons.pt(tightMu[1]));
    double eta1(muons.eta(tightMu[1]));
    skimWeight = 1. - \
      (1. - event_.getScaleFactorMap(Event::Mu23_TrkIsoVVL).get(pt0, eta0)) * \
      (1. - event_.getScaleFactorMap(Event::Mu12_TrkIsoVVL).get(pt0, eta0)) * \
      (1. - event_.getScaleFactorMap(Event::Mu23_TrkIsoVVL).get(pt1, eta1)) * \
      (1. - event_.getScaleFactorMap(Event::Mu12_TrkIsoVVL).get(pt1, eta1));
  }

  return true;
}

fakeskim::DimuonPhotonSkim::DimuonPhotonSkim(Event const& event) :
  FakeSkim(event)
{
  triggers_.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
  triggers_.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
}

fakeskim::DimuonPhotonSkim::DimuonPhotonSkim(DimuonPhotonSkim const& orig) :
  FakeSkim(orig)
{
  triggers_.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
  triggers_.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
}

void
fakeskim::DimuonPhotonSkim::bookBranches(TTree* outTree)
{
  outTree->Branch("tightMu0_2MG", tightMu.data(), "tightMu0_2MG/i");
  outTree->Branch("tightMu1_2MG", tightMu.data() + 1, "tightMu1_2MG/i");
  outTree->Branch("mmm_2MG", &mmm, "mmm_2MG/F");
  outTree->Branch("Photon_2MG_minDRMu", Photon_minDRMu, "Photon_2MG_minDRMu[nPhoton]/F");
  outTree->Branch("Photon_2MG_mmmg", Photon_mmmg, "Photon_2MG_mmmg[nPhoton]/F");
  outTree->Branch("Photon_2MG_njet", Photon_njet, "Photon_2MG_njet[nElectron]/i");
  outTree->Branch("Photon_2MG_nbjet", Photon_nbjet, "Photon_2MG_nbjet[nElectron]/i");
  outTree->Branch("skimWeight_2MG", &skimWeight, "skimWeight_2MG/F");
}

bool
fakeskim::DimuonPhotonSkim::passSkim()
{
  auto& muons(event_.muons);

  unsigned nM(muons.size());
  
  tightMu = {nM, nM};

  for (unsigned iM(0); iM != nM; ++iM) {
    if (!event_.Muon_isTight[iM])
      continue;

    if (tightMu[0] == nM)
      tightMu[0] = iM;
    else {
      if (tightMu[1] != nM || muons.charge(iM) == muons.charge(tightMu[0]))
        return false;
      
      tightMu[1] = iM;
    }
  }

  if (tightMu[1] == nM)
    return false;

  auto pmm(muons.p4(tightMu[0]) + muons.p4(tightMu[1]));
  mmm = pmm.M();
  if (mmm < 20.)
    return false;

  auto& photons(event_.photons);
  unsigned nP(photons.size());

  std::fill_n(Photon_minDRMu, nP, 0.);
  std::fill_n(Photon_mmmg, nP, 0.);
  std::fill_n(Photon_njet, nP, 0);
  std::fill_n(Photon_nbjet, nP, 0);

  bool hasPhoton(false);

  for (unsigned iP(0); iP != nP; ++iP) {
    if ((photons.cutBasedBitmap(iP) & 2) == 0)
      continue;

    hasPhoton = true;

    Photon_minDRMu[iP] = std::sqrt(std::min(photons.deltaR2(iP, muons, tightMu[0]), photons.deltaR2(iP, muons, tightMu[1])));
    Photon_mmmg[iP] = (pmm + photons.p4(iP)).M();
  }

  if (!hasPhoton)
    return false;

  auto& jets(event_.jets);
  unsigned nJ(jets.size());
  unsigned njet(0);
  unsigned nbjet(0);
  unsigned jetIdxMu0(muons.jetIdx(tightMu[0]));
  unsigned jetIdxMu1(muons.jetIdx(tightMu[1]));
  for (unsigned iJ(0); iJ != nJ; ++iJ) {
    if (iJ == jetIdxMu0 || iJ == jetIdxMu1)
      continue;

    if (jets.pt(iJ) < 30. || std::abs(jets.eta(iJ)) > 4.7)
      continue;

    ++njet;

    if (jets.btagDeepB(iJ) > 0.4941)
      ++nbjet;
  }

  for (unsigned iP(0); iP != nP; ++iP) {
    if ((photons.cutBasedBitmap(iP) & 2) == 0)
      continue;

    Photon_njet[iP] = njet;
    Photon_nbjet[iP] = nbjet;

    unsigned jetIdx(photons.jetIdx(iP));
    if (jetIdx < nJ && jetIdx != jetIdxMu0 && jetIdx != jetIdxMu1 && jets.pt(jetIdx) > 30. && std::abs(jets.eta(jetIdx)) < 4.7) {
      --Photon_njet[iP];
      if (jets.btagDeepB(jetIdx) > 0.4941)
        --Photon_nbjet[iP];
    }
  }

  if (event_.isMC()) {
    double pt0(muons.pt(tightMu[0]) * 1.5);
    double eta0(muons.eta(tightMu[0]));
    double pt1(muons.pt(tightMu[1]) * 1.5);
    double eta1(muons.eta(tightMu[1]));
    skimWeight = \
      (1. - (1. - event_.getScaleFactorMap(Event::Mu23_TrkIsoVVL).get(pt0, eta0)) * (1. - event_.getScaleFactorMap(Event::Mu12_TrkIsoVVL).get(pt0, eta0))) * \
      (1. - (1. - event_.getScaleFactorMap(Event::Mu23_TrkIsoVVL).get(pt1, eta1)) * (1. - event_.getScaleFactorMap(Event::Mu12_TrkIsoVVL).get(pt1, eta1)));
  }

  return true;
}

fakeskim::SameSignMuonElectronSkim::SameSignMuonElectronSkim(Event const& event) :
  FakeSkim(event)
{
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  triggers_.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
}

fakeskim::SameSignMuonElectronSkim::SameSignMuonElectronSkim(SameSignMuonElectronSkim const& orig) :
  FakeSkim(orig)
{
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  triggers_.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
}

void
fakeskim::SameSignMuonElectronSkim::bookBranches(TTree* outTree)
{
  outTree->Branch("tightMu_SSME", &tightMu, "tightMu_SSME/i");
  outTree->Branch("skimWeight_SSME", &skimWeight, "skimWeight_SSME/F");
  outTree->Branch("njet_SSME", &njet, "njet_SSME/i");
  outTree->Branch("nbjet_SSME", &nbjet, "nbjet_SSME/i");
}

bool
fakeskim::SameSignMuonElectronSkim::passSkim()
{
  auto& muons(event_.muons);

  unsigned nM(muons.size());
  
  tightMu = nM;

  for (unsigned iM(0); iM != nM; ++iM) {
    if (!event_.Muon_isTight[iM])
      continue;

    if (tightMu != nM)
      return false;
    
    tightMu = iM;
  }

  if (tightMu == nM)
    return false;

  auto& trigObjs(event_.trigObjs);
  std::vector<unsigned> ele23Objs;
  std::vector<unsigned> ele12Objs;
  if (!event_.isMC()) {
    unsigned nO(trigObjs.size());
    for (unsigned iO(0); iO != nO; ++iO) {
      if (trigObjs.id(iO) != 11)
        continue;
      int filterBits(trigObjs.filterBits(iO));
      double pt(trigObjs.pt(iO));
      if ((filterBits & 1) != 0) {
        if (pt > 23.)
          ele23Objs.push_back(iO);
        if (pt > 12.)
          ele12Objs.push_back(iO);
      }
    }
  }

  auto& electrons(event_.electrons);
  unsigned nE(electrons.size());

  int iBaseline(-1);

  for (unsigned iE(0); iE != nE; ++iE) {
    if (!event_.Electron_isBaseline[iE])
      continue;

    if (iBaseline != -1) {
      iBaseline = -1;
      break;
    }

    if (electrons.charge(iE) != muons.charge(tightMu))
      return false;

    iBaseline = iE;

    // Here the requirement is exactly one baseline electron, which has the same charge as the muon, and matches the trigger object.
    if (!event_.isMC()) {
      if (!electrons.deltaR2Match(iE, trigObjs, ele23Objs, 0.1) && !electrons.deltaR2Match(iE, trigObjs, ele12Objs, 0.1))
        return false;
    }
  }

  if (iBaseline == -1)
    return false;

  auto& jets(event_.jets);
  unsigned nJ(jets.size());
  njet = 0;
  nbjet = 0;
  unsigned jetIdxMu(muons.jetIdx(tightMu));
  unsigned jetIdxE(electrons.jetIdx(iBaseline));
  for (unsigned iJ(0); iJ != nJ; ++iJ) {
    if (iJ == jetIdxMu || iJ == jetIdxE)
      continue;

    if (jets.pt(iJ) < 30. || std::abs(jets.eta(iJ)) > 4.7)
      continue;

    ++njet;

    if (jets.btagDeepB(iJ) > 0.4941)
      ++nbjet;
  }

  if (event_.isMC()) {
    double pt(muons.pt(tightMu));
    double eta(muons.eta(tightMu));
    
    skimWeight = 1. - \
      (1. - event_.getScaleFactorMap(Event::Mu23_TrkIsoVVL).get(pt, eta)) * \
      (1. - event_.getScaleFactorMap(Event::Mu12_TrkIsoVVL).get(pt, eta));
  }

  return true;
}

fakeskim::JetElectronSkim::JetElectronSkim(Event const& event) :
  FakeSkim(event)
{
  triggers_.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30");
}

fakeskim::JetElectronSkim::JetElectronSkim(JetElectronSkim const& orig) :
  FakeSkim(orig)
{
  triggers_.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30");
}

void
fakeskim::JetElectronSkim::bookBranches(TTree* outTree)
{
  outTree->Branch("Electron_JE_minDPhiJet", Electron_minDPhiJet, "Electron_JE_minDPhiJet[nElectron]");
  outTree->Branch("Electron_JE_ptOppJet", Electron_ptOppJet, "Electron_JE_ptOppJet[nElectron]");
  outTree->Branch("Electron_JE_njet", Electron_njet, "Electron_JE_njet[nElectron]/i");
  outTree->Branch("Electron_JE_nbjet", Electron_nbjet, "Electron_JE_nbjet[nElectron]/i");
}
  
bool
fakeskim::JetElectronSkim::passSkim()
{
  if (event_.met.pt() > 20.)
    return false;

  auto& trigObjs(event_.trigObjs);
  std::vector<unsigned> ele12Objs;
  if (!event_.isMC()) {
    unsigned nO(trigObjs.size());
    for (unsigned iO(0); iO != nO; ++iO) {
      if (trigObjs.id(iO) != 11)
        continue;
      int filterBits(trigObjs.filterBits(iO));
      double pt(trigObjs.pt(iO));
      if ((filterBits & 1) != 0) {
        if (pt > 12.)
          ele12Objs.push_back(iO);
      }
    }
  }

  auto& electrons(event_.electrons);
  auto& jets(event_.jets);
  unsigned nE(electrons.size());
  unsigned nJ(jets.size());

  std::fill_n(Electron_minDPhiJet, nE, 0.);
  std::fill_n(Electron_ptOppJet, nE, 0.);

  int iElectron(-1);

  for (unsigned iE(0); iE != nE; ++iE) {
    if (!event_.Electron_isBaseline[iE])
      continue;

    if (iElectron != -1) {
      iElectron = -1;
      break;
    }

    if (!event_.isMC()) {
      if (!electrons.deltaR2Match(iE, trigObjs, ele12Objs, 0.1))
        return false;
    }

    iElectron = iE;
  }

  if (iElectron == -1)
    return false;

  unsigned njet(0);
  unsigned nbjet(0);
  double phi(electrons.phi(iElectron));
  unsigned jetIdx(electrons.jetIdx(iElectron));

  double minDPhi(-1.);
  double maxDPhi(-1.);
  double ptMaxDPhi(0.);
  for (unsigned iJ(0); iJ != nJ; ++iJ) {
    if (iJ == jetIdx)
      continue;

    if (jets.pt(iJ) > 30. && std::abs(jets.eta(iJ)) < 4.7) {
      ++njet;

      if (jets.btagDeepB(iJ) > 0.4941)
        ++nbjet;
    }
      
    double dPhi(std::abs(TVector2::Phi_mpi_pi(jets.phi(iJ) - phi)));
    if (minDPhi < 0. || dPhi < minDPhi)
      minDPhi = dPhi;
    if (dPhi > maxDPhi) {
      maxDPhi = dPhi;
      ptMaxDPhi = jets.pt(iJ);
    }
  }

  Electron_minDPhiJet[iElectron] = minDPhi;
  Electron_ptOppJet[iElectron] = ptMaxDPhi;
  Electron_njet[iElectron] = njet;
  Electron_nbjet[iElectron] = nbjet;

  return true;
}

fakeskim::ElectronSkim::ElectronSkim(Event const& event) :
  FakeSkim(event)
{
  triggers_.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  triggers_.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
  triggers_.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30");
}

fakeskim::ElectronSkim::ElectronSkim(ElectronSkim const& orig) :
  FakeSkim(orig)
{
  triggers_.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
  triggers_.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  triggers_.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
  triggers_.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30");
}

void
fakeskim::ElectronSkim::bookBranches(TTree* outTree)
{
  outTree->Branch("Electron_E_isCaloIdLTrackIdLIsoVL", Electron_isCaloIdLTrackIdLIsoVL, "Electron_E_isCaloIdLTrackIdLIsoVL[nElectron]/O");
}

bool
fakeskim::ElectronSkim::passSkim()
{
  auto& trigObjs(event_.trigObjs);
  std::vector<unsigned> ele23Objs;
  std::vector<unsigned> ele12Objs;
  if (!event_.isMC()) {
    unsigned nO(trigObjs.size());
    for (unsigned iO(0); iO != nO; ++iO) {
      if (trigObjs.id(iO) != 11)
        continue;
      int filterBits(trigObjs.filterBits(iO));
      double pt(trigObjs.pt(iO));
      if ((filterBits & 1) != 0) {
        if (pt > 23.)
          ele23Objs.push_back(iO);
        if (pt > 12.)
          ele12Objs.push_back(iO);
      }
    }
  }

  auto& electrons(event_.electrons);
  unsigned nE(electrons.size());

  std::fill_n(Electron_isCaloIdLTrackIdLIsoVL, nE, false);

  bool hasElectron(false);

  for (unsigned iE(0); iE != nE; ++iE) {
    if (!event_.Electron_isBaseline[iE])
      continue;

    hasElectron = true;

    Electron_isCaloIdLTrackIdLIsoVL[iE] = electrons.deltaR2Match(iE, trigObjs, ele23Objs, 0.1) || electrons.deltaR2Match(iE, trigObjs, ele12Objs, 0.1);
  }

  return hasElectron;
}
