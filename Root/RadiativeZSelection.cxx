#include "brendlingerEgammaTestbed/RadiativeZSelection.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"

#include "xAODTruth/xAODTruthHelpers.h"

// this is needed to distribute the algorithm to the workers
ClassImp(RadiativeZSelection)



RadiativeZSelection::RadiativeZSelection(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



RadiativeZSelection::~RadiativeZSelection()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode RadiativeZSelection::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.

  histoStore()->createTH1F("m_lly",      100, 20, 140);
  histoStore()->createTH1F("m_lle",      100, 20, 140);
  histoStore()->createTH1F("m_llegamma", 100, 20, 140);
  histoStore()->createTH1F("m_ll" ,      100, 20, 140);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode RadiativeZSelection::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();

  xAOD::PhotonContainer photons_all = photonHandler()->getCorrectedContainer();
  xAOD::PhotonContainer photons(SG::VIEW_ELEMENTS);
  for (auto ph : photons_all) {
    if (photonHandler()->passPtEtaCuts(ph)) photons.push_back(ph);
  }

  xAOD::ElectronContainer elecs_all = electronHandler()->getCorrectedContainer();
  xAOD::ElectronContainer elecs(SG::VIEW_ELEMENTS);
  for (auto ele : elecs_all) {
    if (electronHandler()->passPtEtaCuts(ele)) elecs.push_back(ele);
  }

  xAOD::MuonContainer      all_muons    = muonHandler()->getCorrectedContainer();
  xAOD::MuonContainer      muons        = muonHandler()->applySelection(all_muons);
  muons.sort(muonHandler()->comparePt);

  int index_mu1 = -1;
  int index_mu2 = -1;
  AssignZbosonIndices(muons,index_mu1,index_mu2);
  if (index_mu1 < 0 || index_mu2 < 0) return EL::StatusCode::SUCCESS;
  if (photons.size() + elecs.size() == 0) return EL::StatusCode::SUCCESS;

  TLorentzVector Zboson = muons[index_mu1]->p4() + muons[index_mu2]->p4();
  TLorentzVector Zboson_withPhoton;
  TLorentzVector candidate_tlv;

  bool candidate_in_electron_container = HighestPtCandidateIsElectron(photons,elecs);
  if (candidate_in_electron_container) {
    Zboson_withPhoton = Zboson + elecs[0]->p4();
    histoStore()->fillTH1F("m_lle", Zboson_withPhoton.M()/HG::GeV);
    candidate_tlv = elecs[0]->p4();
  }
  else {
    Zboson_withPhoton = Zboson + photons[0]->p4();
    histoStore()->fillTH1F("m_lly", Zboson_withPhoton.M()/HG::GeV);
    candidate_tlv = photons[0]->p4();
  }

  if (candidate_in_electron_container) PrintTruthInformation(*elecs[0]);
  else PrintTruthInformation(*photons[0]);

  // std::cout << Zboson.M() << std::endl;
  histoStore()->fillTH1F("m_llegamma", Zboson_withPhoton.M()/HG::GeV);
  histoStore()->fillTH1F("m_ll", Zboson.M()/HG::GeV);

  return EL::StatusCode::SUCCESS;
}

void RadiativeZSelection::AssignZbosonIndices(xAOD::MuonContainer& muons,
                                              int& SFOS_lep1i,
                                              int& SFOS_lep2i){
  double closest_to=91188.;
  double min_delta = 999999999;
  for (unsigned int i=0;i<muons.size();++i) {
    for (unsigned int j=0;j<muons.size();++j) {
      if (j == i) continue;
      if (muons[i]->charge() == muons[j]->charge()) continue;
      if (muons[i]->pt() < muons[j]->pt()) continue;
      TLorentzVector tmp = muons[i]->p4() + muons[j]->p4();
      if (fabs(tmp.M()-closest_to) < min_delta) {
        min_delta = fabs(tmp.M()-closest_to);
        SFOS_lep1i = i;
        SFOS_lep2i = j;
      }
    }
  }
  return;
}

bool RadiativeZSelection::HighestPtCandidateIsElectron(xAOD::PhotonContainer& photons,
                                                       xAOD::ElectronContainer& elecs) {
  // Make sure the photons and electrons are sorted at this point!
  photons.sort(photonHandler()->comparePt);
  elecs.sort(electronHandler()->comparePt);

  if (!elecs.size()) return false;
  if (!photons.size() and elecs.size()) return true;
  if (elecs[0]->pt() > photons[0]->pt()) return true;
  return false;
}

bool RadiativeZSelection::IsTruePhoton(xAOD::Electron elec) {
  int truth_type   = elec.auxdata< int >( "truthType"   );
  int truth_origin = elec.auxdata< int >( "truthOrigin" );

  int firstEgMotherTruthType   = elec.auxdata< int >( "firstEgMotherTruthType"   );
  int firstEgMotherTruthOrigin = elec.auxdata< int >( "firstEgMotherTruthOrigin" );

  // Sherpa Z->eegamma
  if (truth_type == 4  && truth_origin == 3) return true;
  if (firstEgMotherTruthType == 14 && firstEgMotherTruthOrigin == 3) return true;

  // PowhegPythia Z->ee(gamma)
  if (truth_type == 4 && truth_origin == 5 &&
      firstEgMotherTruthType == 15 && firstEgMotherTruthOrigin == 40) return true;

  return false;
}

bool RadiativeZSelection::IsTruePhoton(xAOD::Photon photon) {

  int truth_type   = photon.auxdata< int >( "truthType"   );
  int truth_origin = photon.auxdata< int >( "truthOrigin" );

  // Sherpa Z->eegamma
  if (truth_type == 14 && truth_origin == 3) return true;
  // PowhegPythia Z->ee(gamma)
  if (truth_type == 15 && truth_origin == 40) return true;

  return false;
}

void RadiativeZSelection::PrintTruthInformation(xAOD::Photon& photon) {

  int truth_type,truth_origin,firstEgMotherTruthType,firstEgMotherTruthOrigin;
  int lastEgMotherTruthType,lastEgMotherTruthOrigin;
  int author;

  truth_type = photon.auxdata< int >( "truthType" );
  truth_origin = photon.auxdata< int >( "truthOrigin" );
  firstEgMotherTruthType = photon.auxdata< int >( "firstEgMotherTruthType" );
  firstEgMotherTruthOrigin = photon.auxdata< int >( "firstEgMotherTruthOrigin" );
  lastEgMotherTruthType = photon.auxdata< int >( "lastEgMotherTruthType" );
  lastEgMotherTruthOrigin = photon.auxdata< int >( "lastEgMotherTruthOrigin" );
  author = photon.author();

  // mother pdgID ?
  const xAOD::TruthParticle* eg_candidate_truth = xAOD::TruthHelpers::getTruthParticle(photon);
  int mc_ph_mother_pdgId = -99999;
  int pdgid = -99999;
  if (eg_candidate_truth) {
    pdgid = eg_candidate_truth->pdgId();
    if (eg_candidate_truth->nParents() > 0) {
      const xAOD::TruthParticle* photon_parent = eg_candidate_truth->parent( 0 );
      mc_ph_mother_pdgId = photon_parent->pdgId();
    }
  }

  std::cout << "Candidate: "
            << " photon "
            << " pt: " << photon.pt()
            << " author: " << author
            << " pdgid: "<< pdgid
            << " truthType: " << truth_type
            << " truthOrigin: " << truth_origin
            << " mother PdgID: " << mc_ph_mother_pdgId
            << " firstEgMotherTT: " << firstEgMotherTruthType
            << " firstEgMotherTO: " << firstEgMotherTruthOrigin
            << " lastEgMotherTT: " << lastEgMotherTruthType
            << " lastEgMotherTO: " << lastEgMotherTruthOrigin
            << " Result: " << IsTruePhoton(photon)
            << std::endl;

  return;
}

void RadiativeZSelection::PrintTruthInformation(xAOD::Electron& elec) {

  int truth_type,truth_origin,firstEgMotherTruthType,firstEgMotherTruthOrigin;
  int lastEgMotherTruthType,lastEgMotherTruthOrigin;
  int author;

  truth_type = elec.auxdata< int >( "truthType" );
  truth_origin = elec.auxdata< int >( "truthOrigin" );
  firstEgMotherTruthType = elec.auxdata< int >( "firstEgMotherTruthType" );
  firstEgMotherTruthOrigin = elec.auxdata< int >( "firstEgMotherTruthOrigin" );
  lastEgMotherTruthType = elec.auxdata< int >( "lastEgMotherTruthType" );
  lastEgMotherTruthOrigin = elec.auxdata< int >( "lastEgMotherTruthOrigin" );
  author = elec.author();

  // mother pdgID ?
  const xAOD::TruthParticle* eg_candidate_truth = xAOD::TruthHelpers::getTruthParticle(elec);
  int mc_ph_mother_pdgId = -99999;
  int pdgid = -99999;
  if (eg_candidate_truth) {
    pdgid = eg_candidate_truth->pdgId();
    if (eg_candidate_truth->nParents() > 0) {
      const xAOD::TruthParticle* photon_parent = eg_candidate_truth->parent( 0 );
      mc_ph_mother_pdgId = photon_parent->pdgId();
    }
  }

  std::cout << "Candidate: "
            << " electron "
            << " pt: " << elec.pt()
            << " author: " << author
            << " pdgid: "<< pdgid
            << " truthType: " << truth_type
            << " truthOrigin: " << truth_origin
            << " mother PdgID: " << mc_ph_mother_pdgId
            << " firstEgMotherTT: " << firstEgMotherTruthType
            << " firstEgMotherTO: " << firstEgMotherTruthOrigin
            << " lastEgMotherTT: " << lastEgMotherTruthType
            << " lastEgMotherTO: " << lastEgMotherTruthOrigin
            << " Result: " << IsTruePhoton(elec)
            << std::endl;

  return;
}
