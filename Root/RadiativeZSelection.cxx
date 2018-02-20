#include "brendlingerEgammaTestbed/RadiativeZSelection.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"

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

  xAOD::PhotonContainer photons = photonHandler()->getCorrectedContainer();

  xAOD::ElectronContainer elecs = electronHandler()->getCorrectedContainer();

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
  if (HighestPtCandidateIsElectron(photons,elecs)) {
    Zboson_withPhoton = Zboson + elecs[0]->p4();
    histoStore()->fillTH1F("m_lle", Zboson_withPhoton.M()/HG::GeV);
  }
  else {
    Zboson_withPhoton = Zboson + photons[0]->p4();
    histoStore()->fillTH1F("m_lly", Zboson_withPhoton.M()/HG::GeV);
  }

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
