
#ifndef _MyUnweightedClass_h_
#define _MyUnweightedClass_h_

#include "TRef.h"
#include "TObject.h"
#include "TRefArray.h"
#include "TLorentzVector.h"

#include "classes/SortableObject.h"

class UnweightedParticle: public SortableObject
{
public:
  Int_t PID; // particle HEP ID number | hepevt.idhep[number]

  Int_t Status; // particle status | hepevt.isthep[number]

  Int_t Mother1; 
  Int_t Mother2; 

  Int_t ColorLine1; 
  Int_t ColorLine2; 

  Float_t Px; // particle momentum vector (x component) | hepevt.phep[number][0]
  Float_t Py; // particle momentum vector (y component) | hepevt.phep[number][1]
  Float_t Pz; // particle momentum vector (z component) | hepevt.phep[number][2]
  Float_t E; // particle energy | hepevt.phep[number][3]
  Float_t M; // particle mass

  Float_t PT; // particle transverse momentum
  Float_t Eta; // particle pseudorapidity
  Float_t Phi; // particle azimuthal angle

  Float_t Rapidity; // particle rapidity

  Float_t LifeTime; 
  Float_t Spin; 

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4();

  ClassDef(UnweightedParticle, 1)
};
#endif
