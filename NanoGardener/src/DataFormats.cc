#include "LatinoAnalysis/NanoGardener/interface/DataFormats.h"

namespace nanoaod {
  std::array<BranchConf, ParticleCollection::nUniqueBranches> const ParticleCollection::branchConfs{{
    {'F', "pt"},
    {'F', "eta"},
    {'F', "phi"}
  }};

  std::array<BranchConf, BaseLeptonCollection::nUniqueBranches> const BaseLeptonCollection::branchConfs{{
    {'I', "pdgId"},
    {'F', "dxy"},
    {'F', "dz"},
    {'I', "charge"}
  }};

  std::array<BranchConf, ElectronCollection::nUniqueBranches> const ElectronCollection::branchConfs{{
    {'O', "mvaFall17Iso_WP90"},
    {'I', "cutBased"},
    {'F', "pfRelIso03_all"},
    {'O', "convVeto"},
    {'F', "sieie"},
    {'F', "hoe"},
    {'F', "eInvMinusPInv"},
    {'F', "eCorr"},
    {'F', "deltaEtaSC"},
    {'b', "lostHits"},
    {'I', "jetIdx"},
    {'I', "vidNestedWPBitmap"}
  }};

  std::array<BranchConf, MuonCollection::nUniqueBranches> const MuonCollection::branchConfs{{
    {'O', "tightId"},
    {'O', "isPFcand"},
    {'F', "pfRelIso04_all"},
    {'I', "nStations"},
    {'I', "jetIdx"}
  }};

  std::array<BranchConf, PhotonCollection::nUniqueBranches> const PhotonCollection::branchConfs{{
    {'I', "cutBasedBitmap"},
    {'O', "electronVeto"},
    {'I', "jetIdx"}
  }};

  std::array<BranchConf, LeptonCollection::nUniqueBranches> const LeptonCollection::branchConfs{{
    {'I', "electronIdx"},
    {'I', "muonIdx"},
    {'I', "isLoose"}
  }};

  std::array<BranchConf, JetCollection::nUniqueBranches> const JetCollection::branchConfs{{
    {'I', "jetId"},
    {'I', "puId"},
    {'F', "btagDeepB"},
    {'F', "btagDeepC"},
    {'I', "partonFlavour"},
    {'I', "hadronFlavour"}
  }};

  std::array<BranchConf, TrigObjCollection::nUniqueBranches> const TrigObjCollection::branchConfs{{
    {'I', "id"},
    {'I', "filterBits"}
  }};

  std::array<BranchConf, GenPartCollection::nUniqueBranches> const GenPartCollection::branchConfs{{
    {'F', "mass"},
    {'I', "pdgId"},
    {'I', "status"},
    {'I', "statusFlags"},
    {'I', "genPartIdxMother"}
  }};

  std::array<BranchConf, LHEPartCollection::nUniqueBranches> const LHEPartCollection::branchConfs{{
    {'I', "pdgId"}
  }};

  std::array<BranchConf, MET::nUniqueBranches> const MET::branchConfs{{
    {'F', "pt"},
    {'F', "phi"}
  }};

  void
  Collection::setTreeReader(TTreeReader* reader, TString const& coll)
  {
    reset();

    n_ = std::make_unique<UIntValueReader>(*reader, "n" + coll);
    for (auto& conf : branchConfs_) {
      TString branchName(coll + "_" + conf.second);
      if (reader->GetTree()->GetBranch(branchName) == nullptr)
        readers_.emplace_back(nullptr);
      else {
        switch (conf.first) {
        case 'I':
          readers_.emplace_back(new IntArrayReader(*reader, branchName));
          break;
        case 'F':
          readers_.emplace_back(new FloatArrayReader(*reader, branchName));
          break;
        case 'O':
          readers_.emplace_back(new BoolArrayReader(*reader, branchName));
          break;
        case 'b':
          readers_.emplace_back(new UCharArrayReader(*reader, branchName));
          break;
        default:
          throw std::runtime_error("Unknown branch type!");
        }
      }
    }
  }

  void
  Singleton::setTreeReader(TTreeReader* reader, TString const& obj)
  {
    reset();

    for (auto& conf : branchConfs_) {
      TString branchName(obj + "_" + conf.second);
      if (reader->GetTree()->GetBranch(branchName) == nullptr)
        readers_.emplace_back(nullptr);
      else {
        switch (conf.first) {
        case 'I':
          readers_.emplace_back(new IntValueReader(*reader, branchName));
          break;
        case 'F':
          readers_.emplace_back(new FloatValueReader(*reader, branchName));
          break;
        case 'O':
          readers_.emplace_back(new BoolValueReader(*reader, branchName));
          break;
        case 'b':
          readers_.emplace_back(new UCharValueReader(*reader, branchName));
          break;
        default:
          throw std::runtime_error("Unknown branch type!");
        }
      }
    }
  }

  bool
  ParticleCollection::deltaR2Match(unsigned iP, std::vector<std::pair<double, double>> const& etaPhis, double coneSize) const
  {
    for (auto& etaPhi : etaPhis) {
      if (deltaR2(iP, etaPhi) < coneSize)
        return true;
    }

    return false;
  }

  bool
  ParticleCollection::deltaR2Match(unsigned iP, ParticleCollection const& testColl, double coneSize) const
  {
    for (unsigned iT(0); iT != testColl.size(); ++iT) {
      if (deltaR2(iP, testColl, iT) < coneSize)
        return true;
    }

    return false;
  }

  bool
  ParticleCollection::deltaR2Match(unsigned iP, ParticleCollection const& testColl, std::vector<unsigned> const& iTs, double coneSize) const
  {
    for (unsigned iT : iTs) {
      if (deltaR2(iP, testColl, iT) < coneSize)
        return true;
    }

    return false;
  }

  bool
  ParticleCollection::deltaR2Match(unsigned iP, ParticleCollection const& testColl, bool const* mask, double coneSize) const
  {
    for (unsigned iT(0); iT != testColl.size(); ++iT) {
      if (mask[iT] && deltaR2(iP, testColl, iT) < coneSize)
        return true;
    }

    return false;
  }

  bool
  ParticleCollection::deltaR2PtMatch(unsigned iP, ParticleCollection const& testColl, std::vector<unsigned> const& iTs, double coneSize, double relPtMargin) const
  {
    for (unsigned iT : iTs) {
      if (deltaR2(iP, testColl, iT) > coneSize)
        continue;

      double relPt(pt(iP) / testColl.pt(iT));
      if (relPt > 1. - relPtMargin && relPt < 1. + relPtMargin)
        return true;
    }

    return false;
  }

  bool
  ElectronCollection::passCuts(unsigned _iE, WorkingPoint _wp, unsigned _cuts) const
  {
    unsigned bitMap(vidNestedWPBitmap(_iE));

    for (unsigned iC(0); iC != nCuts; ++iC) {
      if (((_cuts >> iC) & 0x1) != 0 && ((bitMap >> (iC * 3)) & 0x7) < _wp)
        return false;
    }

    return true;
  }

  void
  LeptonCollection::setTreeReader(TTreeReader* reader, TString const& coll)
  {
    Collection::setTreeReader(reader, coll);
    for (auto& te : isTightElectron_)
      te.second = std::make_unique<IntArrayReader>(*reader, coll + "_isTightElectron_" + te.first);
    for (auto& tm : isTightMuon_)
      tm.second = std::make_unique<IntArrayReader>(*reader, coll + "_isTightMuon_" + tm.first);
  }

  void
  LeptonCollection::reset()
  {
    for (auto& te : isTightElectron_)
      te.second.reset();
    for (auto& tm : isTightMuon_)
      tm.second.reset();
    ParticleCollection::reset();
  }

  bool
  LeptonCollection::isTight(unsigned iP) const
  {
    if (std::abs(pdgId(iP)) == 11) {
      for (auto& te : isTightElectron_) {
        if (te.second->At(iP) != 0)
          return true;
      }
    }
    else {
      for (auto& tm : isTightMuon_) {
        if (tm.second->At(iP) != 0)
          return true;
      }
    }
    return false;
  }

}
