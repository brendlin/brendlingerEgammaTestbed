#ifndef brendlingerEgammaTestbed_RadiativeZSelection_H
#define brendlingerEgammaTestbed_RadiativeZSelection_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

class RadiativeZSelection : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  // Tree *myTree; //!
  // TH1 *myHist; //!



public:
  // this is a standard constructor
  RadiativeZSelection() { }
  RadiativeZSelection(const char *name);
  virtual ~RadiativeZSelection();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();


  // this is needed to distribute the algorithm to the workers
  ClassDef(RadiativeZSelection, 1);
};

#endif // brendlingerEgammaTestbed_RadiativeZSelection_H