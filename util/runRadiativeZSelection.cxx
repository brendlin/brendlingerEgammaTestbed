#include "brendlingerEgammaTestbed/RadiativeZSelection.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  RadiativeZSelection *alg = new RadiativeZSelection("RadiativeZSelection");

  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
