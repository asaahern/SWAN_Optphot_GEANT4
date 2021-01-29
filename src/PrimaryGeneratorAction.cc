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
//`s
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "CLHEP/Units/PhysicalConstants.h"


using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0)
{
    // define particle type and energy
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");

    G4double wavelength = G4RandGauss::shoot(170.*nm, 6.*nm); // Xe scintillation
    G4double energy = h_Planck * c_light / wavelength;

    // create particle gun
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // source parameters
    G4double source_position    = 21.45*mm; //20.5*mm;  --> 2019 data, should be = cylinder_length - 2*Cu_thickness
    G4double source_radius      = 5.*mm; //12.5*mm;  --> 2019 data
    G4double pressure           = 10.;
    G4double alpha_range        = 24./pressure;

    // define initial position:
    // isotropic random distribution within a semi-sphere of radius = alpha range
    G4double theta  = std::acos(G4UniformRand())*rad;
    G4double phi    = G4UniformRand() * 360.0*deg;
    G4double rho    = G4UniformRand() * alpha_range;
    G4double posX   = rho * std::sin(theta) * std::cos(phi);
    G4double posY   = rho * std::sin(theta) * std::sin(phi);
    G4double posZ   = rho * std::cos(theta);
    // centered following uniform random distribution from circumference containing source
    G4double theta_s= G4UniformRand() * 360.0*deg;
    G4double rho_s  = std::sqrt(G4UniformRand() * source_radius * source_radius);
    G4double posX_s = rho_s * std::cos(theta_s);
    G4double posY_s = rho_s * std::sin(theta_s);
    fParticleGun->SetParticlePosition(G4ThreeVector(posX_s+posX, posY_s+posY, source_position-posZ));
    //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, 21.4));  // launch all photons from point source //

    // define momentum direction: isotropic and random
    theta           = std::acos(1-2* G4UniformRand())*rad;
    phi             = G4UniformRand() * 360.0*deg;
    G4double mx     = std::sin(theta) * std::cos(phi);
    G4double my     = std::sin(theta) * std::sin(phi);
    G4double mz     = std::cos(theta);
    G4ThreeVector momentumDirection = G4ThreeVector(mx, my, mz);
    fParticleGun->SetParticleMomentumDirection(momentumDirection);
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1)); // launch all photons directed towards base

    // define polarization direction: random but contained in the plane perpendicular to the momentum direction of each photon
    G4ThreeVector normal(1., 0., 0.);
    G4ThreeVector product = normal.cross(momentumDirection);
    G4double modul2 = product * product;

    G4ThreeVector e_perpend(0., 0., 1.);
    if (modul2 > 0.) e_perpend = (1. / std::sqrt(modul2)) * product;
    G4ThreeVector e_paralle = e_perpend.cross(momentumDirection);

    G4double angle = G4UniformRand() * 360.0 * deg;
    G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
    fParticleGun->SetParticlePolarization(polar);

}

