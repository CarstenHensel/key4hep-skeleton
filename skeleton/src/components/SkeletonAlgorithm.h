/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once

// GAUDI
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "k4FWCore/DataHandle.h"
#include "podio/UserDataCollection.h"

#include "TFile.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <string>
#include <vector>

class SkeletonAlgorithm : public Gaudi::Algorithm {
public:
  explicit SkeletonAlgorithm(const std::string&, ISvcLocator*);
  virtual ~SkeletonAlgorithm();
  /**  Initialize.
   *   @return status code
   */
  StatusCode initialize() override;
  /**  Execute.
   *   @return status code
   */
  StatusCode execute(const EventContext&) const override;
  /**  Finalize.
   *   @return status code
   */
  StatusCode finalize() override;

  // write to tree
  void fillEvent() const; // call once per event

  // Helper method to add an MC particle
  void addMCParticle(float pt, float eta, float phi, float e, int pdgId, int status, int motherPdgId) const;
  void addMCParticle(const TLorentzVector& p4, int pdgId, int status, int motherPdgId) const;

  // Public variables to be filled before calling fillEvent()
  mutable int runNumber;
  mutable int eventNumber;
  mutable float lumiWeight;
  mutable float weight;
  mutable float sqs_eff;
  // TODO not sure what to do with this one
  // I'll keep it for historically reasons
  mutable int channelType; // 1 = leptonic Z, 2 = hadronic Z
  mutable int proc_id;
  mutable int nElectrons;
  mutable int nMuons;
  mutable int nJets;
  mutable bool passElectronPreSelection;
  mutable bool passMuonPreSelection;
  mutable bool passJetPreSelection;

  // leptons
  // electrons
  mutable float ele_pt[2];
  mutable float ele_eta[2];
  mutable float ele_phi[2];
  mutable float ele_e[2];
  mutable int ele_charge[2];
  mutable float di_electron_mass;
  mutable float di_electron_pt;
  mutable float ele_delta_phi;
  mutable float ele_delta_r;

  // muons
  mutable float muon_pt[2];
  mutable float muon_eta[2];
  mutable float muon_phi[2];
  mutable float muon_e[2];
  mutable int muon_charge[2];
  mutable float di_muon_mass;
  mutable float di_muon_pt;
  mutable float muon_delta_phi;
  mutable float muon_delta_r;

  // jets
  mutable float jet_pt[8];
  mutable float jet_eta[8];
  mutable float jet_phi[8];
  mutable float jet_e[8];
  mutable float jet_m[8];

  // MET/Recoil
  mutable float MET;
  mutable float MET_phi;
  mutable float MET_px;
  mutable float MET_py;
  mutable float recoil_mass;
  mutable float recoil_pt;
  mutable float recoil_cosTheta;

  // kinematics
  mutable float acoplanarity;
  mutable float acollinearity;
  mutable float deltaPhi;
  mutable float deltaR;
  mutable float visible_mass;
  mutable float visible_energy;
  mutable float visible_pt;

  // MC Truth
  mutable int nMCParticles;
  mutable std::vector<float> mc_pt;
  mutable std::vector<float> mc_eta;
  mutable std::vector<float> mc_phi;
  mutable std::vector<float> mc_e;
  mutable std::vector<int> mc_pdgId;
  mutable std::vector<int> mc_status;
  mutable std::vector<int> mc_motherPdgId;

private:
  // member variable
  mutable k4FWCore::DataHandle<edm4hep::ReconstructedParticleCollection> m_recoParticleCollHandle{
      "ReconstructedParticleCollection", Gaudi::DataHandle::Reader, this};

  mutable k4FWCore::DataHandle<edm4hep::ReconstructedParticleCollection> m_isolatedLeptonsCollHandle{
      "IsolatedLeptonsCollection", Gaudi::DataHandle::Reader, this};

  mutable k4FWCore::DataHandle<edm4hep::EventHeaderCollection> m_eventHeaderCollHandle{"EventHeaderCollection",
                                                                             Gaudi::DataHandle::Reader, this};

  mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_mcParticleCollHandle{"MCParticleColl", Gaudi::DataHandle::Reader,
                                                                           this};

  mutable k4FWCore::DataHandle<edm4hep::ReconstructedParticleCollection> m_jetFinderCollHandle{
      "JetFinderColl", // property name (what you set in options)
      Gaudi::DataHandle::Reader, this};

  mutable k4FWCore::DataHandle<podio::UserDataCollection<float>> m_outMET{"MET", Gaudi::DataHandle::Writer, this};

  // addiing some code to extract information from the event header
  mutable k4FWCore::DataHandle<edm4hep::EventHeaderCollection> m_eventHeader{"EventHeader", Gaudi::DataHandle::Reader, this};

  Gaudi::Property<std::vector<std::string>> m_outputs{this, "Outputs", {}, "Output collection names"};

  mutable std::string m_selection;
  mutable double luminosity_weight;
  mutable int processID;
  mutable int n_events_generated;
  mutable float targetLumi;
  mutable float cross_section;
  mutable std::string processName;
  mutable std::string root_output_file;

  mutable int m_event_counter = 0;

  mutable int m_member = 0;

  // tree related
  TFile* outFile;
  TTree* tree;

  void setupBranches();

  // outsourcing MET calculation to a method
  // Usage: auto [MET, MET_px, MET_py, MET_phi] = this->getMET(...)
  std::tuple<double, double, double, double> getMET(const edm4hep::ReconstructedParticleCollection* particles) const;

  // outsourcing the calculation of visible kinematics
  struct VisibleKinematics {
    double mass;   // invariant mass of visible particles
    double pt;     // transverse momentum
    double energy; // total energy
  };

  VisibleKinematics computeVisibleKinematics(const edm4hep::ReconstructedParticleCollection* particles) const;

  // fill electron and muon variables
  void fillLeptons(const edm4hep::ReconstructedParticleCollection* particles) const;

  // fill jet brances
  void fillJets(const edm4hep::ReconstructedParticleCollection* jets) const;

  // some helper methods to calculate dela_phi, delta_r, and di_particle_pt
  double getDeltaPhi(double phi1, double phi2) const;
  double getDeltaR(const edm4hep::ReconstructedParticle& p1, const edm4hep::ReconstructedParticle& p2) const;
  double diParticlePt(const edm4hep::ReconstructedParticle& p1, const edm4hep::ReconstructedParticle& p2) const;
};
