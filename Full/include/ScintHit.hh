//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: LXeScintHit.hh 72250 2013-07-12 08:59:26Z gcosmo $
//
/// \file optical/LXe/include/LXeScintHit.hh
/// \brief Definition of the LXeScintHit class
//
//
#ifndef ScintHit_h
#define ScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

#include "tls.hh"

class ScintHit : public G4VHit
{
  public:
 
    ScintHit();
    ScintHit(G4VPhysicalVolume* pVol);
    virtual ~ScintHit();
    ScintHit(const ScintHit &right);
    const ScintHit& operator=(const ScintHit &right);
    G4int operator==(const ScintHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
 
    virtual void Draw();
    virtual void Print();

    inline void SetEdep(G4double de) { fEdep = de; }
    inline void AddEdep(G4double de) { fEdep += de; }
    inline G4double GetEdep() { return fEdep; }

    inline void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    inline G4ThreeVector GetPos() { return fPos; }

    inline const G4VPhysicalVolume * GetPhysV() { return fPhysVol; }

  private:
    G4double fEdep;
    G4ThreeVector fPos;
    const G4VPhysicalVolume* fPhysVol;

};

typedef G4THitsCollection<ScintHit> ScintHitsCollection;

extern G4ThreadLocal G4Allocator<ScintHit>* ScintHitAllocator;

inline void* ScintHit::operator new(size_t)
{
  if(!ScintHitAllocator)
      ScintHitAllocator = new G4Allocator<ScintHit>;
  return (void *) ScintHitAllocator->MallocSingle();
}

inline void ScintHit::operator delete(void *aHit)
{
  ScintHitAllocator->FreeSingle((ScintHit*) aHit);
}

#endif