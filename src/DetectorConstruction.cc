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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
	: G4VUserDetectorConstruction()
{	// World
	world_half_side = 100*mm;

	// TPC chamber
	// Source support
	al_radius       = 12.5*mm;
	// Copper anode
	cu_radius       = 32.5*mm;
	cu_thickness    = 1.*mm;

	// Teflon cylinder
	cylinder_length = 23.45*mm; //5.*mm;

	// Drift to base
	drift_length    = 26.5*mm;

	// Multi Wire Proportional Chambers
	MW_in_radius    = 30.*mm;
	MW_out_radius   = 49.*mm;
	MW_thickness    = 3.75*mm;

	// TPC base
	// Teflon base
	base_radius    = 50.*mm;
	base_thickness = 2.5*mm;

	// Photomultipliers
	PMT_radius     = 12.5*mm;
	PMT_thickness  = 2.5*mm;
	phCa_radius    = 12.5*mm;
	phc_thickness  = 0.05*mm;
	jacket_length  = 5.*mm;

	disable_drift_wall = false;
	disable_al_support = false;
	disable_cu_anode   = false;
	disable_PMT_base   = false;
	disable_PMT_jackets = false;
	disable_MWPC       = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() { ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = true;

	G4NistManager* nist = G4NistManager::Instance();

	
	// ------------ Generate & Add Material Properties Table ------------


	// Xe ------------------------------------------

	// gas density
	G4double pressure = 10.;
	G4double temperature = 298.15 * kelvin;
	G4Material* matXe = nist->ConstructNewGasMaterial("GXe", "G4_Xe", temperature, pressure);
	G4double density_xe = (pressure / atmosphere) * 131.29 / (temperature / kelvin * 82.058); // g/cm^3

	// refractive index calculations
	const G4int nXe_entries = 100;
	G4double nXe_energy[nXe_entries];
	G4double refractiveIndex_Xe[nXe_entries];
	G4double energy_start = 1.8 * eV;
	G4double energy_finish = 8.9 * eV;

	for (G4int i = 0; i < nXe_entries; i++) {
		nXe_energy[i] = (energy_start + i * (energy_finish - energy_start) / nXe_entries);


		// Formula for the refractive index taken from
		// A. Baldini et al., "Liquid Xe scintillation calorimetry 
		// and Xe optical properties", arXiv:physics/0401072v1 [physics.ins-det]

		// The Lorentz-Lorenz equation (also known as Clausius-Mossotti equation)
		// relates the refractive index of a fluid with its density:
		// (n^2 - 1) / (n^2 + 2) = - A · d_M,     (1)
		// where n is the refractive index, d_M is the molar density and
		// A is the first refractivity viral coefficient:
		// A(E) = \sum_i^3 P_i / (E^2 - E_i^2),   (2)
		// with:
		G4double P[3] = { 71.23, 77.75, 1384.89 }; // [eV^3 cm3 / mole]
		G4double E[3] = { 8.4, 8.81, 13.2 };       // [eV]

		// Note.- Equation (1) has, actually, a sign difference with respect 
		// to the one appearing in the reference. Otherwise, it yields values
		// for the refractive index below 1.

		// Let's calculate the virial coefficient.
		// We won't use the implicit system of units of Geant4 because
		// it results in loss of numerical precision.

		G4double energy_ite = nXe_energy[i] / eV;

		G4double virial = 0.;
		for (G4int j = 0; j < 3; j++)
			virial = virial + P[j] / (energy_ite * energy_ite - E[j] * E[j]);

		G4double mol_density = density_xe / 131.29;
		G4double alpha = virial * mol_density;

		// Isolating now the n2 from equation (1) and taking the square root
		refractiveIndex_Xe[i] = (1. - 2 * alpha) / (1. + alpha);

		if (refractiveIndex_Xe[i] < 1.) {
			// "Non-physical refractive index for energy "
			refractiveIndex_Xe[i] = 1.;
		}

	}
	assert(sizeof(refractiveIndex_Xe) == sizeof(nXe_energy));

	G4MaterialPropertiesTable* Xe_MPT1 = new G4MaterialPropertiesTable();

	Xe_MPT1->AddProperty("RINDEX", nXe_energy, refractiveIndex_Xe, nXe_entries)
		->SetSpline(true);

	G4cout << "Xe G4MaterialPropertiesTable" << G4endl;

	matXe->SetMaterialPropertiesTable(Xe_MPT1);


	// PMT window ----------------------------------------
	G4Material* PMT_window_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

	// define window refractive index (source: janis.com/Libraries/Window_Transmissions/FusedSilicaUVGrade_SiO2_TransmissionCurveDataSheet.sflb.ashx)
	G4double photonEnergy[] = { 7.3074*eV, 6.7149*eV, 6.2113*eV, 5.7941*eV, 4.4319*eV, 4.1121*eV, 3.4035*eV,
								3.0704*eV, 2.8505*eV, 2.2748*eV, 2.1141*eV, 2.108*eV, 1.9296*eV, 1.8928*eV, 1.441*eV };
	
	G4double refractiveIndex1[] = { 1.6150, 1.5750, 1.5500, 1.5337, 1.4940, 1.4872, 1.4745,
									1.4696, 1.4666, 1.4601, 1.4585, 1.4584, 1.4567, 1.4564, 1.4525 };

	const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

	G4MaterialPropertiesTable* pmtMPT1 = new G4MaterialPropertiesTable();

	pmtMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries);
	PMT_window_mat->SetMaterialPropertiesTable(pmtMPT1);


	// PMT photocathode --------------------------------
	G4Material* PMT_phc_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");


	// Copper anode  -----------------------------------
	G4Material* copper_mat = nist->FindOrBuildMaterial("G4_Cu");


	// Americium source Al layer  ----------------------
	G4Material* source_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");


	// drift wall   ------------------------------------
	G4Material* teflon_mat = nist->FindOrBuildMaterial("G4_TEFLON");

	const G4int num3 = 4;
	// source: C. Silva et al. “Reflectance of Polytetrafluoroethylene (PTFE) for Xenon Scintillation Light”,
	// J. Appl. Phys. 107 (2010) 064902
	G4double ephoton_teflon[num3] = { 2.21*eV , 3.95*eV, 4.87*eV, 7.3*eV };
	//G4double reflectivity_teflon[num3] = { 0.98, 0.93, 0.85, 0.61};
	G4double reflectivity_teflon[num3] = { 0., 0., 0., 0.};


  // ------------ Colours ---------------

	G4VisAttributes * blue = new G4VisAttributes(G4Colour(0. ,0.8 ,0.9, 0.6));
	blue -> SetVisibility(true);
	blue -> SetForceSolid(true);

	G4VisAttributes * black = new G4VisAttributes(G4Colour(0.2 ,0.2 ,0.2));
	black -> SetVisibility(true);
	black -> SetForceSolid(true);

	G4VisAttributes * redCu = new G4VisAttributes(G4Colour(0.3 ,0.1 ,0.));
	redCu -> SetVisibility(true);
	redCu -> SetForceSolid(true);

	G4VisAttributes * graySS = new G4VisAttributes(G4Colour(0.7 ,0.7 ,0.7));
	graySS -> SetVisibility(true);
	graySS -> SetForceSolid(true);

	G4VisAttributes * yellow = new G4VisAttributes(G4Colour(1. ,1. ,0., 1.0));
	yellow -> SetVisibility(true);
	yellow -> SetForceSolid(true);

	G4VisAttributes * orange = new G4VisAttributes(G4Colour(0.6 ,0.4 ,0., 1.0));
	orange -> SetVisibility(true);
	orange -> SetForceSolid(true);

  // ------------- Volumes --------------

    // The experimental Hall
	G4Box* expHall_box = new G4Box("World", world_half_side, world_half_side, world_half_side);
	G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box, matXe, "World", 0, 0, 0);
	G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0, G4ThreeVector(), expHall_log, "World", 0, false, 0, checkOverlaps);

	// PMT windows
	G4Tubs* PMT_win1 = new G4Tubs("PMT_win1", 0.* cm, PMT_radius, PMT_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtWindow1_log = new G4LogicalVolume(PMT_win1, PMT_window_mat, "PMT_win1", 0, 0, 0);
	G4VPhysicalVolume* pmtWindow1_phys = new G4PVPlacement(0, G4ThreeVector(-25.*mm, 0, -cylinder_length - drift_length + phc_thickness*2 + PMT_thickness), pmtWindow1_log, "PMT_win1", expHall_log, false, 0, checkOverlaps);

	G4Tubs* PMT_win2 = new G4Tubs("PMT_win2", 0.* cm, PMT_radius, PMT_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtWindow2_log = new G4LogicalVolume(PMT_win2, PMT_window_mat, "PMT_win2", 0, 0, 0);
	G4VPhysicalVolume* pmtWindow2_phys = new G4PVPlacement(0, G4ThreeVector(25.*mm, 0, -cylinder_length - drift_length + phc_thickness*2 + PMT_thickness), pmtWindow2_log, "PMT_win2", expHall_log, false, 0, checkOverlaps);

	G4Tubs* PMT_win3 = new G4Tubs("PMT_win3", 0.* cm, PMT_radius, PMT_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtWindow3_log = new G4LogicalVolume(PMT_win3, PMT_window_mat, "PMT_win3", 0, 0, 0);
	G4VPhysicalVolume* pmtWindow3_phys = new G4PVPlacement(0, G4ThreeVector(0, -25.*mm, -cylinder_length - drift_length + phc_thickness*2 + PMT_thickness), pmtWindow3_log, "PMT_win3", expHall_log, false, 0, checkOverlaps);

	G4Tubs* PMT_win4 = new G4Tubs("PMT_win4", 0.* cm, PMT_radius, PMT_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtWindow4_log = new G4LogicalVolume(PMT_win4, PMT_window_mat, "PMT_win4", 0, 0, 0);
	G4VPhysicalVolume* pmtWindow4_phys = new G4PVPlacement(0, G4ThreeVector(0, 25.*mm, -cylinder_length - drift_length + phc_thickness*2 + PMT_thickness), pmtWindow4_log, "PMT_win4", expHall_log, false, 0, checkOverlaps);

	pmtWindow1_log -> SetVisAttributes(blue);
	pmtWindow2_log -> SetVisAttributes(blue);
	pmtWindow3_log -> SetVisAttributes(blue);
	pmtWindow4_log -> SetVisAttributes(blue);
	

	// PMT photocathode
	G4Tubs* PMT_phc1 = new G4Tubs("PMT_phCa1", 0. * cm, phCa_radius, phc_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtPhc1_log = new G4LogicalVolume(PMT_phc1, PMT_phc_mat, "PMT_phCa1", 0, 0, 0);
	G4VPhysicalVolume* pmtPhc1_phys = new G4PVPlacement(0, G4ThreeVector(0, 25.*mm, -cylinder_length - drift_length + phc_thickness), pmtPhc1_log, "PMT_phCa1", expHall_log, false, 0, checkOverlaps);

	G4Tubs* PMT_phc2 = new G4Tubs("PMT_phCa2", 0. * cm, phCa_radius, phc_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtPhc2_log = new G4LogicalVolume(PMT_phc2, PMT_phc_mat, "PMT_phCa2", 0, 0, 0);
	G4VPhysicalVolume* pmtPhc2_phys = new G4PVPlacement(0, G4ThreeVector(0, -25.*mm, -cylinder_length - drift_length + phc_thickness), pmtPhc2_log, "PMT_phCa2", expHall_log, false, 0, checkOverlaps);

	G4Tubs* PMT_phc3 = new G4Tubs("PMT_phCa3", 0. * cm, phCa_radius, phc_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtPhc3_log = new G4LogicalVolume(PMT_phc3, PMT_phc_mat, "PMT_phCa3", 0, 0, 0);
	G4VPhysicalVolume* pmtPhc3_phys = new G4PVPlacement(0, G4ThreeVector(25.*mm, 0, -cylinder_length - drift_length + phc_thickness), pmtPhc3_log, "PMT_phCa3", expHall_log, false, 0, checkOverlaps);

	G4Tubs* PMT_phc4 = new G4Tubs("PMT_phCa4", 0. * cm, phCa_radius, phc_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtPhc4_log = new G4LogicalVolume(PMT_phc4, PMT_phc_mat, "PMT_phCa4", 0, 0, 0);
	G4VPhysicalVolume* pmtPhc4_phys = new G4PVPlacement(0, G4ThreeVector(-25.*mm, 0, -cylinder_length - drift_length + phc_thickness), pmtPhc4_log, "PMT_phCa4", expHall_log, false, 0, checkOverlaps);

	pmtPhc1_log -> SetVisAttributes(black);
	pmtPhc2_log -> SetVisAttributes(black);
	pmtPhc3_log -> SetVisAttributes(black);
	pmtPhc4_log -> SetVisAttributes(black);

	// Copper anode
	if (!disable_cu_anode) {

		// volume
		G4Tubs* Ca = new G4Tubs("cathode", al_radius, cu_radius, cu_thickness, 0. * deg, 360. * deg);
		G4LogicalVolume* Ca_log = new G4LogicalVolume(Ca, copper_mat, "cathode", 0, 0, 0);
		G4VPhysicalVolume* Ca_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, cylinder_length-cu_thickness), Ca_log, "cathode", expHall_log, false, 0, checkOverlaps);

		// optical surface
		G4OpticalSurface* opCa = new G4OpticalSurface("Ca_Surface");
		G4LogicalSkinSurface* Ca_Surface = new G4LogicalSkinSurface("Ca_Surface", Ca_log, opCa);

		opCa->SetType(dielectric_metal);
		opCa->SetFinish(ground);
		opCa->SetModel(glisur);
		opCa->SetPolish(0.4);

		const G4int num = 11;
		// source: https://www.photonics.com/Articles/Mirrors_Coating_Choice_Makes_a_Difference/a25501
		G4double ephoton_Ca[num] = { 0.12*eV, 0.31*eV, 0.62*eV, 1.24*eV, 1.77*eV, 2.067*eV, 2.48*eV, 3.1*eV, 4.13*eV, 6.2*eV,  8.9*eV };
		G4double reflectivity_Ca[num] = { 0.982, 0.979, 0.97, 0.97, 0.955, 0.85, 0.592, 0.51, 0.346, 0.358, 0.358 };
		//G4double reflectivity_Ca[num] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

		G4MaterialPropertiesTable* CaST2 = new G4MaterialPropertiesTable();
		CaST2->AddProperty("REFLECTIVITY", ephoton_Ca, reflectivity_Ca, num);
		opCa->SetMaterialPropertiesTable(CaST2);

		Ca_log -> SetVisAttributes(redCu);
	}


	// Stainless-steel layer over americium source 
	if (!disable_al_support) {

		// volume
		G4Tubs* deWi = new G4Tubs("detector_window", 0.0, al_radius, cu_thickness, 0. * deg, 360. * deg);
		G4LogicalVolume* deWi_log = new G4LogicalVolume(deWi, source_mat, "detector_window", 0, 0, 0);
		G4VPhysicalVolume* deWi_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, cylinder_length-cu_thickness), deWi_log, "cathode", expHall_log, false, 0, checkOverlaps);

		// optical surface (aluminium evaporated ~70 nm) 
		G4OpticalSurface* op_dWind = new G4OpticalSurface("dWin_Surface");
		G4LogicalSkinSurface* dWind_Surface = new G4LogicalSkinSurface("dWin_Surface", deWi_log, op_dWind);

		op_dWind->SetType(dielectric_metal);
		op_dWind->SetFinish(ground);
		op_dWind->SetModel(glisur);
		op_dWind->SetPolish(0.4);

		const G4int num4 = 19;
		// source: Noble-gas liquid detectors: measurement of light diffusion and reflectivity
		// on commonly adopted inner surface materials , S. Bricola et al (only for Xe 2nd continuum)
		// source: STAINLESS STEEL SOLAR MIRRORS - A MATERIAL FEASIBILITY STUDY , Ame ROOS , fig6 (for the remaining wavelength spectrum)    
		G4double ephoton_win[num4] = {0.35*eV, 0.48*eV, 0.64*eV, 0.79*eV, 0.93*eV, 1.06*eV, 1.22*eV, 1.4*eV, 1.58*eV, 1.81*eV,
									  2.04*eV, 2.27*eV, 2.6*eV, 2.91*eV, 3.21*eV, 3.57*eV, 3.91*eV, 4.13*eV, 7.3*eV};
		G4double reflectivity_win[num4] = {0.854, 0.82, 0.783, 0.765, 0.743, 0.729, 0.712, 0.703, 0.696, 0.683,
										   0.676, 0.659, 0.627, 0.605, 0.576, 0.55, 0.541, 0.541, 0.57};        
		//G4double reflectivity_win[num4] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};


		G4MaterialPropertiesTable* dWindST2 = new G4MaterialPropertiesTable();
		dWindST2->AddProperty("REFLECTIVITY", ephoton_win, reflectivity_win, num4);

		op_dWind->SetMaterialPropertiesTable(dWindST2);

		deWi_log -> SetVisAttributes(graySS);
	}


	// Teflon cylinder
	if (!disable_drift_wall) {

		// volume
		G4Tubs* drift_wall = new G4Tubs("drift_wall", cu_radius, cu_radius + cu_thickness, cylinder_length, 0. * deg, 360. * deg);
		G4LogicalVolume* drift_wall_log	= new G4LogicalVolume(drift_wall, teflon_mat, "drift_wall", 0, 0, 0);
		G4VPhysicalVolume* drift_wall_phys = new G4PVPlacement(0, G4ThreeVector(), drift_wall_log, "drift_wall", expHall_log, false, 0, checkOverlaps);

		// optical surface
		G4OpticalSurface* op_dWall = new G4OpticalSurface("dWall_Surface");
		G4LogicalSkinSurface* dWall_Surface = new G4LogicalSkinSurface("dWall_Surface", drift_wall_log, op_dWall);

		G4MaterialPropertiesTable* dWallST2 = new G4MaterialPropertiesTable();
		dWallST2->AddProperty("REFLECTIVITY", ephoton_teflon, reflectivity_teflon, num3);

		op_dWall->SetMaterialPropertiesTable(dWallST2);
		op_dWall->SetType(dielectric_metal);
		op_dWall->SetFinish(ground);
		op_dWall->SetPolish(0.2);
		op_dWall->SetModel(glisur);

		drift_wall_log  -> SetVisAttributes(yellow);
	}

	// MWPC
	if (!disable_MWPC) {

		// volume
		G4Tubs* MWPC = new G4Tubs("MWPC", MW_in_radius, MW_out_radius, MW_thickness, 0. * deg, 360. * deg);
		G4LogicalVolume* MWPC_log	= new G4LogicalVolume(MWPC, teflon_mat, "MWPC", 0, 0, 0);
		G4VPhysicalVolume* MWPC_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -cylinder_length - MW_thickness), MWPC_log, "MWPC", expHall_log, false, 0, checkOverlaps);

		MWPC_log  -> SetVisAttributes(orange);
	}



	// PMT jackets
	if (!disable_PMT_jackets) {

		// volume
		G4Tubs* PMT_jack1 = new G4Tubs("PMT_jack1", PMT_radius, PMT_radius + cu_thickness, jacket_length, 0. * deg, 360. * deg);
		G4LogicalVolume* PMT_jack1_log	= new G4LogicalVolume(PMT_jack1, teflon_mat, "PMT_jack1", 0, 0, 0);
		G4VPhysicalVolume* PMT_jack1_phys = new G4PVPlacement(0, G4ThreeVector(0, 25.*mm, -cylinder_length - drift_length + jacket_length), PMT_jack1_log, "PMT_jack1", expHall_log, false, 0, checkOverlaps);

		G4Tubs* PMT_jack2 = new G4Tubs("PMT_jack2", PMT_radius, PMT_radius + cu_thickness, jacket_length, 0. * deg, 360. * deg);
		G4LogicalVolume* PMT_jack2_log	= new G4LogicalVolume(PMT_jack2, teflon_mat, "PMT_jack2", 0, 0, 0);
		G4VPhysicalVolume* PMT_jack2_phys = new G4PVPlacement(0, G4ThreeVector(0, -25.*mm, -cylinder_length - drift_length + jacket_length), PMT_jack2_log, "PMT_jack2", expHall_log, false, 0, checkOverlaps);

		//G4Tubs* PMT_jack3 = new G4Tubs("PMT_jack3", PMT_radius, PMT_radius + cu_thickness, jacket_length, 0. * deg, 360. * deg);
		//G4LogicalVolume* PMT_jack3_log	= new G4LogicalVolume(PMT_jack3, teflon_mat, "PMT_jack3", 0, 0, 0);
		//G4VPhysicalVolume* PMT_jack3_phys = new G4PVPlacement(0, G4ThreeVector(25.*mm, 0, -cylinder_length - drift_length + jacket_length), PMT_jack3_log, "PMT_jack3", expHall_log, false, 0, checkOverlaps);

		//G4Tubs* PMT_jack4 = new G4Tubs("PMT_jack4", PMT_radius, PMT_radius + cu_thickness, jacket_length, 0. * deg, 360. * deg);
		//G4LogicalVolume* PMT_jack4_log	= new G4LogicalVolume(PMT_jack4, teflon_mat, "PMT_jack4", 0, 0, 0);
		//G4VPhysicalVolume* PMT_jack4_phys = new G4PVPlacement(0, G4ThreeVector(-25.*mm, 0, -cylinder_length - drift_length + jacket_length), PMT_jack4_log, "PMT_jack4", expHall_log, false, 0, checkOverlaps);

		// optical surface
		G4MaterialPropertiesTable* PMT_jack_R = new G4MaterialPropertiesTable();
		PMT_jack_R->AddProperty("REFLECTIVITY", ephoton_teflon, reflectivity_teflon, num3);

		G4OpticalSurface* op_PMT_jack = new G4OpticalSurface("PMT_jack_Surface");
		op_PMT_jack->SetMaterialPropertiesTable(PMT_jack_R);
		op_PMT_jack->SetType(dielectric_metal);
		op_PMT_jack->SetFinish(ground);
		op_PMT_jack->SetPolish(0.2);
		op_PMT_jack->SetModel(glisur);

		G4LogicalSkinSurface* PMT_jack_Surface1 = new G4LogicalSkinSurface("PMT_jack1_Surface", PMT_jack1_log, op_PMT_jack);
		G4LogicalSkinSurface* PMT_jack_Surface2 = new G4LogicalSkinSurface("PMT_jack2_Surface", PMT_jack2_log, op_PMT_jack);

		PMT_jack1_log  -> SetVisAttributes(yellow);
		PMT_jack2_log  -> SetVisAttributes(yellow);
		//PMT_jack3_log  -> SetVisAttributes(yellow);
		//PMT_jack4_log  -> SetVisAttributes(yellow);
	}

	// PMT base
	if (!disable_PMT_base) {

		// volume
		G4Tubs* PMT_base = new G4Tubs("PMT_base", 0., base_radius, base_thickness, 0. * deg, 360. * deg);
		G4LogicalVolume* PMT_base_log	= new G4LogicalVolume(PMT_base, teflon_mat, "PMT_base", 0, 0, 0);
		G4VPhysicalVolume* PMT_base_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -cylinder_length - drift_length - base_thickness), PMT_base_log, "PMT_base", expHall_log, false, 0, checkOverlaps);

		// optical surface
		G4OpticalSurface* op_PMT_base = new G4OpticalSurface("PMT_base_Surface");
		G4LogicalSkinSurface* PMT_base_Surface = new G4LogicalSkinSurface("PMT_base_Surface", PMT_base_log, op_PMT_base);

		G4MaterialPropertiesTable* PMT_base_ST2 = new G4MaterialPropertiesTable();
		PMT_base_ST2->AddProperty("REFLECTIVITY", ephoton_teflon, reflectivity_teflon, num3);

		op_PMT_base->SetMaterialPropertiesTable(PMT_base_ST2);
		op_PMT_base->SetType(dielectric_metal);
		op_PMT_base->SetFinish(ground);
		op_PMT_base->SetPolish(0.2);
		op_PMT_base->SetModel(glisur);

		PMT_base_log  -> SetVisAttributes(yellow);
	}



	// ------------- Surfaces --------------

	// PMT Quartz Window
	G4OpticalSurface* opPMT_window = new G4OpticalSurface("PMT_Surface");
	opPMT_window->SetType(dielectric_LUTDAVIS);
	opPMT_window->SetFinish(Polished_LUT);
	opPMT_window->SetModel(DAVIS);

	G4LogicalBorderSurface* PMT_Surface1 = new G4LogicalBorderSurface("PMT_BorderSurface1", pmtWindow1_phys, expHall_phys, opPMT_window);
	G4OpticalSurface* PMTopticalSurface1 = dynamic_cast <G4OpticalSurface*>
		(PMT_Surface1->GetSurface(pmtWindow1_phys, expHall_phys)->
			GetSurfaceProperty());
	if (PMTopticalSurface1) PMTopticalSurface1->DumpInfo();

	G4LogicalBorderSurface* PMT_Surface2 = new G4LogicalBorderSurface("PMT_BorderSurface2", pmtWindow2_phys, expHall_phys, opPMT_window);
	G4OpticalSurface* PMTopticalSurface2 = dynamic_cast <G4OpticalSurface*>
		(PMT_Surface1->GetSurface(pmtWindow2_phys, expHall_phys)->
			GetSurfaceProperty());
	if (PMTopticalSurface2) PMTopticalSurface2->DumpInfo();

	G4LogicalBorderSurface* PMT_Surface3 = new G4LogicalBorderSurface("PMT_BorderSurface3", pmtWindow3_phys, expHall_phys, opPMT_window);
	G4OpticalSurface* PMTopticalSurface3 = dynamic_cast <G4OpticalSurface*>
		(PMT_Surface1->GetSurface(pmtWindow3_phys, expHall_phys)->
			GetSurfaceProperty());
	if (PMTopticalSurface3) PMTopticalSurface3->DumpInfo();

	G4LogicalBorderSurface* PMT_Surface4 = new G4LogicalBorderSurface("PMT_BorderSurface4", pmtWindow4_phys, expHall_phys, opPMT_window);
	G4OpticalSurface* PMTopticalSurface4 = dynamic_cast <G4OpticalSurface*>
		(PMT_Surface1->GetSurface(pmtWindow4_phys, expHall_phys)->
			GetSurfaceProperty());
	if (PMTopticalSurface4) PMTopticalSurface4->DumpInfo();

	// photocathode

	const G4int num2 = 2;
	G4double ephoton_abs[num2] = { 1.0 * eV , 10.0 * eV };
	G4double efficiency[num2] = { 1.0 , 1.0 };
	G4double reflectivity_ph[num2] = { 0.0, 0.0 };

	G4OpticalSurface* opPhCa = new G4OpticalSurface("PhCa_Surface");

	G4MaterialPropertiesTable* PhCaST2 = new G4MaterialPropertiesTable();
	PhCaST2->AddProperty("EFFICIENCY", ephoton_abs, efficiency, num2);
	PhCaST2->AddProperty("REFLECTIVITY", ephoton_abs, reflectivity_ph, num2);
	
	opPhCa->SetMaterialPropertiesTable(PhCaST2);
	opPhCa->SetType(dielectric_metal);
	opPhCa->SetFinish(polished);
	opPhCa->SetModel(glisur);

	G4LogicalSkinSurface* PhCa_Surface1 = new G4LogicalSkinSurface("PhCa_SkinSurface1", pmtPhc1_log, opPhCa);
	G4LogicalSkinSurface* PhCa_Surface2 = new G4LogicalSkinSurface("PhCa_SkinSurface2", pmtPhc2_log, opPhCa);
	G4LogicalSkinSurface* PhCa_Surface3 = new G4LogicalSkinSurface("PhCa_SkinSurface3", pmtPhc3_log, opPhCa);
	G4LogicalSkinSurface* PhCa_Surface4 = new G4LogicalSkinSurface("PhCa_SkinSurface4", pmtPhc4_log, opPhCa);
	
	//always return the physical World
	return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
