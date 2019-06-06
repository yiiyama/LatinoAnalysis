#ifndef NANOAOD_DATAFORMATS_h
#define NANOAOD_DATAFORMATS_h

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TString.h"

typedef TTreeReaderArray<Float_t> FloatArrayReader;
typedef TTreeReaderValue<Float_t> FloatValueReader;
typedef TTreeReaderArray<Int_t> IntArrayReader;
typedef TTreeReaderValue<Int_t> IntValueReader;
typedef TTreeReaderArray<Bool_t> BoolArrayReader;
typedef TTreeReaderValue<Bool_t> BoolValueReader;
typedef TTreeReaderArray<UChar_t> UCharArrayReader;
typedef TTreeReaderValue<UChar_t> UCharValueReader;
typedef TTreeReaderArray<UInt_t> UIntArrayReader;
typedef TTreeReaderValue<UInt_t> UIntValueReader;
typedef TTreeReaderValue<ULong64_t> ULong64ValueReader;
typedef std::unique_ptr<FloatArrayReader> FloatArrayReaderPtr;
typedef std::unique_ptr<FloatValueReader> FloatValueReaderPtr;
typedef std::unique_ptr<IntArrayReader> IntArrayReaderPtr;
typedef std::unique_ptr<IntValueReader> IntValueReaderPtr;
typedef std::unique_ptr<BoolArrayReader> BoolArrayReaderPtr;
typedef std::unique_ptr<BoolValueReader> BoolValueReaderPtr;
typedef std::unique_ptr<UCharArrayReader> UCharArrayReaderPtr;
typedef std::unique_ptr<UCharValueReader> UCharValueReaderPtr;
typedef std::unique_ptr<UIntArrayReader> UIntArrayReaderPtr;
typedef std::unique_ptr<UIntValueReader> UIntValueReaderPtr;
typedef std::unique_ptr<ULong64ValueReader> ULong64ValueReaderPtr;

namespace nanoaod {

  typedef std::pair<char, TString> BranchConf;

  class Collection {
  public:
    Collection() {}
    virtual ~Collection() {}
    virtual void setTreeReader(TTreeReader*, TString const&);
    virtual void reset() {
      n_.reset();
      readers_.clear();
    }
    bool valid() const { return bool(n_); }
    unsigned size() const { return *n_->Get(); }

    enum Branch {
      nBranches
    };

    double getF(unsigned br, unsigned iP) const { return static_cast<FloatArrayReader&>(*readers_[br]).At(iP); }
    int getI(unsigned br, unsigned iP) const { return static_cast<IntArrayReader&>(*readers_[br]).At(iP); }
    bool getO(unsigned br, unsigned iP) const { return static_cast<BoolArrayReader&>(*readers_[br]).At(iP); }
    unsigned getb(unsigned br, unsigned iP) const { return static_cast<UCharArrayReader&>(*readers_[br]).At(iP); }

  protected:
    std::vector<BranchConf> branchConfs_{};

  private:
    UIntValueReaderPtr n_{};
    std::vector<std::unique_ptr<ROOT::Internal::TTreeReaderArrayBase>> readers_{};
  };

  class Singleton {
  public:
    Singleton() {}
    virtual ~Singleton() {}
    virtual void setTreeReader(TTreeReader*, TString const&);
    virtual void reset() {
      readers_.clear();
    }

    enum Branch {
      nBranches
    };

    double getF(unsigned br) const { return *static_cast<FloatValueReader&>(*readers_[br]).Get(); }
    int getI(unsigned br) const { return *static_cast<IntValueReader&>(*readers_[br]).Get(); }
    bool getO(unsigned br) const { return *static_cast<BoolValueReader&>(*readers_[br]).Get(); }
    unsigned getb(unsigned br) const { return *static_cast<UCharValueReader&>(*readers_[br]).Get(); }

  protected:
    std::vector<BranchConf> branchConfs_{};

  private:
    std::vector<std::unique_ptr<ROOT::Internal::TTreeReaderValueBase>> readers_{};
  };
  
  class ParticleCollection : public Collection {
  public:
    ParticleCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _pt = Collection::nBranches,
      _eta,
      _phi,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - Collection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    double pt(unsigned iP) const { return getF(_pt, iP); };
    double eta(unsigned iP) const { return getF(_eta, iP); };
    double phi(unsigned iP) const { return getF(_phi, iP); };

    virtual double mass(unsigned iP) const { return 0.; };
    TLorentzVector p4(unsigned iP) const {
      TLorentzVector p;
      p.SetPtEtaPhiM(pt(iP), eta(iP), phi(iP), mass(iP));
      return p;
    }
    TVector2 vpt(unsigned iP) const { TVector2 ret; ret.SetMagPhi(pt(iP), phi(iP)); return ret; }
    double deltaR2(unsigned iP, std::pair<double, double> const& etaPhi) const {
      double dEta(eta(iP) - etaPhi.first);
      double dPhi(TVector2::Phi_mpi_pi(phi(iP) - etaPhi.second));
      return dEta * dEta + dPhi * dPhi;
    }
    double deltaR2(unsigned iP, ParticleCollection const& testColl, unsigned iT) const { return deltaR2(iP, {testColl.eta(iT), testColl.phi(iT)}); }
    bool deltaR2Match(unsigned iP, std::vector<std::pair<double, double>> const&, double coneSize) const;
    bool deltaR2Match(unsigned iP, ParticleCollection const&, double coneSize) const;
    bool deltaR2Match(unsigned iP, ParticleCollection const&, std::vector<unsigned> const& iTs, double coneSize) const;
    bool deltaR2Match(unsigned iP, ParticleCollection const&, bool const* mask, double coneSize) const;
    bool deltaR2PtMatch(unsigned iP, ParticleCollection const&, std::vector<unsigned> const& iTs, double coneSize, double relPtMargin) const;
  };

  class BaseLeptonCollection : public ParticleCollection {
  public:
    BaseLeptonCollection() : ParticleCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _pdgId = ParticleCollection::nBranches,
      _dxy,
      _dz,
      _charge,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - ParticleCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    int pdgId(unsigned iP) const { return getI(_pdgId, iP); }
    double dxy(unsigned iP) const { return getF(_dxy, iP); }
    double dz(unsigned iP) const { return getF(_dz, iP); }
    int charge(unsigned iP) const { return getI(_charge, iP); }
  };

  class ElectronCollection : public BaseLeptonCollection {
  public:
    ElectronCollection() : BaseLeptonCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _mvaFall17Iso_WP90 = BaseLeptonCollection::nBranches,
      _cutBased,
      _pfRelIso03_all,
      _convVeto,
      _sieie,
      _hoe,
      _eInvMinusPInv,
      _eCorr,
      _deltaEtaSC,
      _lostHits,
      _jetIdx,
      _vidNestedWPBitmap,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - BaseLeptonCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    bool mvaFall17Iso_WP90(unsigned iE) const { return getO(_mvaFall17Iso_WP90, iE); }
    int cutBased(unsigned iE) const { return getI(_cutBased, iE); }
    double pfRelIso03_all(unsigned iE) const { return getF(_pfRelIso03_all, iE); }
    bool convVeto(unsigned iE) const { return getO(_convVeto, iE); }
    double sieie(unsigned iE) const { return getF(_sieie, iE); }
    double hoe(unsigned iE) const { return getF(_hoe, iE); }
    double eInvMinusPInv(unsigned iE) const { return getF(_eInvMinusPInv, iE); }
    double eCorr(unsigned iE) const { return getF(_eCorr, iE); }
    double deltaEtaSC(unsigned iE) const { return getF(_deltaEtaSC, iE); }
    unsigned char lostHits(unsigned iE) const { return getb(_lostHits, iE); }
    int jetIdx(unsigned iE) const { return getI(_jetIdx, iE); }
    int vidNestedWPBitmap(unsigned iE) const { return getI(_vidNestedWPBitmap, iE); }

    enum Cut {
      MinPtCut,
      GsfEleSCEtaMultiRangeCut,
      GsfEleDEtaInSeedCut,
      GsfEleDPhiInCut,
      GsfEleFull5x5SigmaIEtaIEtaCut,
      GsfEleHadronicOverEMEnergyScaledCut,
      GsfEleEInverseMinusPInverseCut,
      GsfEleRelPFIsoScaledCut,
      GsfEleConversionVetoCut,
      GsfEleMissingHitsCut,
      nCuts
    };

    enum WorkingPoint {
      Veto = 1,
      Loose,
      Medium,
      Tight,
      nWorkingPoints
    };
  
    bool passCuts(unsigned iE, WorkingPoint wp, unsigned cuts) const;
    double etaSC(unsigned iP) const { return eta(iP) + deltaEtaSC(iP); }
  };

  class MuonCollection : public BaseLeptonCollection {
  public:
    MuonCollection() : BaseLeptonCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _tightId = BaseLeptonCollection::nBranches,
      _isPFcand,
      _pfRelIso04_all,
      _nStations,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - BaseLeptonCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    bool tightId(unsigned iM) const { return getO(_tightId, iM); }
    bool isPFcand(unsigned iM) const { return getO(_isPFcand, iM); }
    double pfRelIso04_all(unsigned iM) const { return getF(_pfRelIso04_all, iM); }
    int nStations(unsigned iM) const { return getI(_nStations, iM); }
  };

  class PhotonCollection : public ParticleCollection {
  public:
    PhotonCollection() : ParticleCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _cutBasedBitmap = ParticleCollection::nBranches,
      _electronVeto,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - ParticleCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    int cutBasedBitmap(unsigned iP) const { return getI(_cutBasedBitmap, iP); }
    bool electronVeto(unsigned iP) const { return getO(_electronVeto, iP); }
  };

  class LeptonCollection : public BaseLeptonCollection {
  public:
    LeptonCollection() : BaseLeptonCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }
    void setTreeReader(TTreeReader*, TString const&) override;
    void reset() override;

    enum Branch {
      _electronIdx = BaseLeptonCollection::nBranches,
      _muonIdx,
      _isLoose,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - BaseLeptonCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    int electronIdx(unsigned iP) const { return getI(_electronIdx, iP); };
    int muonIdx(unsigned iP) const { return getI(_muonIdx, iP); };
    int isLoose(unsigned iP) const { return getI(_isLoose, iP); }
    bool isTight(unsigned iP) const;

    void addTightElectronWP(TString const& wp) { isTightElectron_.emplace(wp, nullptr); }
    void addTightMuonWP(TString const& wp) { isTightMuon_.emplace(wp, nullptr); }

  private:
    std::map<TString, IntArrayReaderPtr> isTightElectron_{};
    std::map<TString, IntArrayReaderPtr> isTightMuon_{};
  };

  class JetCollection : public ParticleCollection {
  public:
    JetCollection() : ParticleCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _jetId = ParticleCollection::nBranches,
      _puId,
      _btagDeepB,
      _btagDeepC,
      _partonFlavour,
      _hadronFlavour,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - ParticleCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    int jetId(unsigned iJ) const { return getI(_jetId, iJ); }
    int puId(unsigned iJ) const { return getI(_puId, iJ); }
    double btagDeepB(unsigned iJ) const { return getF(_btagDeepB, iJ); }
    double btagDeepC(unsigned iJ) const { return getF(_btagDeepC, iJ); }
    int partonFlavour(unsigned iJ) const { return getI(_partonFlavour, iJ); }
    int hadronFlavour(unsigned iJ) const { return getI(_hadronFlavour, iJ); }
  };

  class TrigObjCollection : public ParticleCollection {
  public:
    TrigObjCollection() : ParticleCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _id = ParticleCollection::nBranches,
      _filterBits,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - ParticleCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    int id(unsigned iO) const { return getI(_id, iO); };
    int filterBits(unsigned iO) const { return getI(_filterBits, iO); };

    bool hasFilterBits(unsigned iO, unsigned bits) const { return (filterBits(iO) & bits) != 0; }
  };

  class GenPartCollection : public ParticleCollection {
  public:
    GenPartCollection() : ParticleCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _mass = ParticleCollection::nBranches,
      _pdgId,
      _status,
      _statusFlags,
      _genPartIdxMother,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - ParticleCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    double mass(unsigned iP) const override { return getF(_mass, iP); };
    int pdgId(unsigned iP) const { return getI(_pdgId, iP); };
    int status(unsigned iP) const { return getI(_status, iP); }
    int statusFlags(unsigned iP) const { return getI(_statusFlags, iP); }
    int genPartIdxMother(unsigned iP) const { return getI(_genPartIdxMother, iP); }

    enum StatusBit {
      kIsPrompt = 0,
      kIsDecayedLeptonHadron,
      kIsTauDecayProduct,
      kIsPromptTauDecayProduct,
      kIsDirectTauDecayProduct,
      kIsDirectPromptTauDecayProduct,
      kIsDirectHadronDecayProduct,
      kIsHardProcess,
      kFromHardProcess,
      kIsHardProcessTauDecayProduct,
      kIsDirectHardProcessTauDecayProduct,
      kFromHardProcessBeforeFSR,
      kIsFirstCopy,
      kIsLastCopy,
      kIsLastCopyBeforeFSR,
      nStatusBits
    };

    bool checkStatus(unsigned iP, StatusBit b) { return (statusFlags(iP) & (1 << b)) != 0; }
    bool checkStatus(unsigned iP, StatusBit b1, StatusBit b2) { return (statusFlags(iP) & ((1 << b1) | (1 << b2))) != 0; }
  };

  class LHEPartCollection : public ParticleCollection {
  public:
    LHEPartCollection() : ParticleCollection() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _pdgId = ParticleCollection::nBranches,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - ParticleCollection::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    int pdgId(unsigned iP) const { return getI(_pdgId, iP); };
  };

  class MET : public Singleton {
  public:
    MET() : Singleton() {
      branchConfs_.insert(branchConfs_.end(), branchConfs.begin(), branchConfs.end());
    }

    enum Branch {
      _pt = Singleton::nBranches,
      _phi,
      nBranches
    };
    static constexpr unsigned nUniqueBranches{nBranches - Singleton::nBranches};
    static std::array<BranchConf, nUniqueBranches> const branchConfs;

    double pt() const { return getF(_pt); }
    double phi() const { return getF(_phi); }

    TVector2 vpt() const { TVector2 ret; ret.SetMagPhi(pt(), phi()); return ret; }
  };

}

#endif
