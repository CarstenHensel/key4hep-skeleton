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
#include <TLorentzVector.h>

#include <vector>

#include "GaudiKernel/MsgStream.h"
#include "SelectEvents.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/Vector3d.h"
#include "podio/Frame.h"
#include "podio/UserDataCollection.h"

DECLARE_COMPONENT(SelectEvents)

SelectEvents::SelectEvents(const std::string& aName, ISvcLocator* aSvcLoc) : Gaudi::Algorithm(aName, aSvcLoc) {
  declareProperty("RecoParticleColl", m_recoParticleCollHandle, "RecoParticle collection");
  declareProperty("IsolatedLeptonsColl", m_isolatedLeptonsCollHandle, "Isolated Leptons collection");
  declareProperty("EventHeaderColl", m_eventHeaderCollHandle, "Event Header collection");
  declareProperty("MCParticleColl", m_mcParticleCollHandle, "MC Particle collection");
  declareProperty("JetFinderColl", m_jetFinderCollHandle, "Jet Finder Collection");

  declareProperty("Outputs", m_outMET, "Name of the output MET collection");

  declareProperty("cross_section", cross_section, "cross_section");
  declareProperty("processID", processID, "processID");
  declareProperty("n_events_generated", n_events_generated, "n_events_generated");
  declareProperty("targetLumi", targetLumi, "targetLumi");
  declareProperty("processName", processName, "processName");
  declareProperty("root_output_file", root_output_file, "root_output_file");
}

SelectEvents::~SelectEvents() {}

StatusCode SelectEvents::initialize() {
  m_event_counter = 0;

  // set up output file
  std::string filename = root_output_file;
  outFile = new TFile(filename.c_str(), "RECREATE");

  // creat output tree
  tree = new TTree("events", "Higgs to Invisible Analysis Tree");
  tree->SetDirectory(0);  // This prevents ROOT from automatically deleting it
  this->setupBranches();

  // set some global variables
  lumiWeight = cross_section * targetLumi / n_events_generated;

  // set selection
  m_selection = "mumu"; // TODO: promote variable to property
  return StatusCode::SUCCESS;
}

StatusCode SelectEvents::execute(const EventContext& event) const {

  // increasing event counter
  m_event_counter += 1;

  // getting the pointers to the collections I need
  const auto* isoLeptonColl = m_isolatedLeptonsCollHandle.get();
  const auto* recoColl = m_recoParticleCollHandle.get();
  const auto* eventHeaderColl = m_eventHeaderCollHandle.get();
  const auto* mcParticleColl = m_mcParticleCollHandle.get();

  // testing the jetfinder coll:

  const auto* jetFinderColl = m_jetFinderCollHandle.get();
  if (!jetFinderColl) {
    error() << "No jet collection found at " << m_jetFinderCollHandle << endmsg;
    return StatusCode::FAILURE;
  }

  // this part is not needed anymore if I keep it to one algorithm only
  // keeping it just in case
  auto userMET = m_outMET.createAndPut();
  userMET->push_back(50.0 + m_event_counter);

  // read and set event variables
  if (eventHeaderColl && !eventHeaderColl->empty()) {
    eventNumber = eventHeaderColl->at(0).getEventNumber();
    runNumber = eventHeaderColl->at(0).getRunNumber();
  }
  sqs_eff = 250.0;   // TODO: find a way to calculate effective center-of-mass energy
  proc_id = processID;

  // get MET
  auto [met, met_px, met_py, met_phi] = this->getMET(recoColl);
  MET = met;
  MET_px = met_px;
  MET_py = met_py;
  MET_phi = met_phi;

  // compute visible kinematics
  VisibleKinematics vis = computeVisibleKinematics(recoColl);
  visible_mass = vis.mass;
  visible_pt = vis.pt;
  visible_energy = vis.energy;

  // fill lepton kinematics
  this->fillLeptons(recoColl);

  // fill jet kinematics
  this->fillJets(jetFinderColl);

  // splitting up selection
  if (m_selection == "ee") {
    // TODO: set up ee selection
  } else if (m_selection == "mumu") {
    // TODO: set up mumu selection
  } else if (m_selection == "jj") {
    // TODO: set up jetjet selection
  }

  // MC particles
  for (const auto& mc : *mcParticleColl) {
    auto parents = mc.getParents();
    int parentPDG;
    if (!parents.empty()) {
      const auto& parent = parents[0]; // if you just want the first parent
      parentPDG = parent.getPDG();
    } else {
      parentPDG = 0;
    }

    TLorentzVector p4;
    auto mom = mc.getMomentum();
    p4.SetPxPyPzE(mom.x, mom.y, mom.z, mc.getEnergy());
    addMCParticle(p4, mc.getPDG(), mc.getSimulatorStatus(), parentPDG);
  }

  fillEvent();
  return StatusCode::SUCCESS;
}

StatusCode SelectEvents::finalize() {
  if (outFile && tree) {
    outFile->cd();
    tree->Write();
    outFile->Close();
  }

  info() << "Before delete - tree pointer: " << tree << endmsg;
  if (tree) {
      info() << "Tree name: " << tree->GetName() << endmsg;
      info() << "Tree entries: " << tree->GetEntries() << endmsg;
  }
    
  if (tree) {
      delete tree;
      tree = nullptr;
      info() << "Tree deleted successfully" << endmsg;
  } else {
      info() << "Tree was already null" << endmsg;
  }

  delete outFile;
  outFile = nullptr;
  
  info() << "SelectEvents algorithm finished!" << endmsg;
  return StatusCode::SUCCESS;
}

// calculate MET
std::tuple<double, double, double, double>
SelectEvents::getMET(const edm4hep::ReconstructedParticleCollection* particles) const {
  // Usage: auto [MET, MET_px, MET_py, MET_phi] = this->getMET(...)
  // calculate MET
  float sum_px = 0.0f;
  float sum_py = 0.0f;

  for (const auto& p : *particles) {
    sum_px += p.getMomentum().x;
    sum_py += p.getMomentum().y;
  }

  float met_px = -sum_px;
  float met_py = -sum_py;
  float met = std::sqrt(met_px * met_px + met_py * met_py);
  float met_phi = std::atan2(met_py, met_px);

  return {met, met_px, met_py, met_phi};
}

// Function to compute visible kinematics
SelectEvents::VisibleKinematics
SelectEvents::computeVisibleKinematics(const edm4hep::ReconstructedParticleCollection* particles) const {
  TLorentzVector total(0, 0, 0, 0);

  for (const auto& p : *particles) {
    const auto& mom = p.getMomentum();
    TLorentzVector v(mom.x, mom.y, mom.z, p.getEnergy());
    total += v;
  }

  VisibleKinematics result;
  result.mass = total.M();
  result.pt = total.Pt();
  result.energy = total.E();
  return result;
}

// fill lepton variable
void SelectEvents::fillLeptons(const edm4hep::ReconstructedParticleCollection* particles) const {

  std::vector<edm4hep::ReconstructedParticle> electrons;
  std::vector<edm4hep::ReconstructedParticle> muons;

  // Loop over reconstructed particles
  for (const auto& p : *particles) {
    // Simple loose preselection (can tune)
    float pt = hypot(p.getMomentum()[0], p.getMomentum()[1]);
    // TODO: double check loose pre-selction
    if (pt > 5.0 && fabs(p.getMomentum()[2] / p.getEnergy()) < 0.98) {
      if (abs(p.getType()) == 11) {
        electrons.push_back(p);
      } else if (abs(p.getType()) == 13) {
        muons.push_back(p);
      }
    }
  }

  // Sort by descending energy
  std::sort(electrons.begin(), electrons.end(), [](auto& a, auto& b) { return a.getEnergy() > b.getEnergy(); });
  std::sort(muons.begin(), muons.end(), [](auto& a, auto& b) { return a.getEnergy() > b.getEnergy(); });

  // set tree variablaes
  nElectrons = electrons.size();
  for (int i = 0; i < 2; ++i) {
    if (i < nElectrons) {
      ele_pt[i] = hypot(electrons[i].getMomentum()[0], electrons[i].getMomentum()[1]);
      ele_eta[i] = atanh(electrons[i].getMomentum()[2] / electrons[i].getEnergy());
      ele_phi[i] = atan2(electrons[i].getMomentum()[1], electrons[i].getMomentum()[0]);
      ele_e[i] = electrons[i].getEnergy();
      ele_charge[i] = electrons[i].getCharge();
    } else {
      ele_pt[i] = ele_eta[i] = ele_phi[i] = ele_e[i] = -99;
      ele_charge[i] = 0;
    }
  }

  nMuons = muons.size();
  for (int i = 0; i < 2; ++i) {
    if (i < nMuons) {
      muon_pt[i] = hypot(muons[i].getMomentum()[0], muons[i].getMomentum()[1]);
      muon_eta[i] = atanh(muons[i].getMomentum()[2] / muons[i].getEnergy());
      muon_phi[i] = atan2(muons[i].getMomentum()[1], muons[i].getMomentum()[0]);
      muon_e[i] = muons[i].getEnergy();
      muon_charge[i] = muons[i].getCharge();
    } else {
      muon_pt[i] = muon_eta[i] = muon_phi[i] = muon_e[i] = -99;
      muon_charge[i] = 0;
    }
  }

  if (nElectrons >= 2) {
    // the '+' operator is not defined in edm4hep::Vector3f
    // so I need a work around, not elegant, but so what
    const auto& p1 = electrons[0].getMomentum();
    const auto& p2 = electrons[1].getMomentum();
    double e1 = electrons[0].getEnergy();
    double e2 = electrons[1].getEnergy();

    TLorentzVector p4;
    p4.SetPxPyPzE(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z, e1 + e2);

    di_electron_mass = p4.M();
    ele_delta_phi = this->getDeltaPhi(ele_phi[0], ele_phi[1]);
    ele_delta_r = this->getDeltaR(electrons[0], electrons[1]);
  } else {
    di_electron_mass = -99.0;
    ele_delta_phi = -99.0;
    ele_delta_r = -99.0;
  }

  if (nMuons >= 2) {
    // the '+' operator is not defined in edm4hep::Vector3f
    // so I need a work around, not elegant, but so what
    const auto& p1 = muons[0].getMomentum();
    const auto& p2 = muons[1].getMomentum();
    double e1 = muons[0].getEnergy();
    double e2 = muons[1].getEnergy();

    TLorentzVector p4;
    p4.SetPxPyPzE(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z, e1 + e2);

    di_muon_mass = p4.M();
    muon_delta_phi = this->getDeltaPhi(muon_phi[0], muon_phi[1]);
    muon_delta_r = this->getDeltaR(muons[0], muons[1]);
  } else {
    di_muon_mass = -99.0;
    muon_delta_phi = -99.0;
    muon_delta_r = -99.0;
  }
}

// fill the jet branches
void SelectEvents::fillJets(const edm4hep::ReconstructedParticleCollection* jetColl) const {

  std::vector<edm4hep::ReconstructedParticle> jets;
  for (const auto& j : *jetColl) {
    jets.push_back(j);
  }
  // sort jets according to their energy
  std::sort(jets.begin(), jets.end(), [](auto& a, auto& b) { return a.getEnergy() > b.getEnergy(); });

  nJets = jets.size();
  // Loop over jets
  for (int i = 0; i < nJets; i++) {

    // Four-momentum
    auto p4 = jets[i].getMomentum(); // returns edm4hep::Vector3f
    double energy = jets[i].getEnergy();
    double px = p4.x;
    double py = p4.y;
    double pz = p4.z;
    double pt = std::sqrt(px * px + py * py);
    double mass = jets[i].getMass();

    // Compute η, φ
    double phi = std::atan2(py, px);
    double p = std::sqrt(px * px + py * py + pz * pz);
    double eta = 0.5 * std::log((p + pz) / (p - pz + 1e-9)); // protect div by 0

    // Fill vectors
    jet_pt[i] = pt;
    jet_e[i] = energy;
    jet_m[i] = mass;
    jet_eta[i] = eta;
    jet_phi[i] = phi;

    // Jet constituents (optional, depending on if you need them)
    // jet_nConstituents.push_back(jet.getParticles().size());
  }
}

// helper methods
// wrap deltaPhi to [-pi, pi]
double SelectEvents::getDeltaPhi(double phi1, double phi2) const {
  double dphi = phi1 - phi2;
  while (dphi > M_PI)
    dphi -= 2 * M_PI;
  while (dphi <= -M_PI)
    dphi += 2 * M_PI;
  return dphi;
}

double SelectEvents::getDeltaR(const edm4hep::ReconstructedParticle& p1, const edm4hep::ReconstructedParticle& p2) const {
  // compute pseudorapidities

  // another missing functionality:
  // I need to calculate the vector magnitude by hand

  const auto& mom1 = p1.getMomentum();
  double mag1 = std::sqrt(mom1.x * mom1.x + mom1.y * mom1.y + mom1.z * mom1.z);
  const auto& mom2 = p2.getMomentum();
  double mag2 = std::sqrt(mom2.x * mom2.x + mom2.y * mom2.y + mom2.z * mom2.z);
  double theta1 = std::acos(p1.getMomentum().z / mag1);
  double theta2 = std::acos(p2.getMomentum().z / mag2);

  double eta1 = -std::log(std::tan(theta1 / 2.0));
  double eta2 = -std::log(std::tan(theta2 / 2.0));

  // compute delta phi
  double dphi = getDeltaPhi(std::atan2(p1.getMomentum().y, p1.getMomentum().x),
                            std::atan2(p2.getMomentum().y, p2.getMomentum().x));

  return std::sqrt((eta1 - eta2) * (eta1 - eta2) + dphi * dphi);
}

double SelectEvents::diParticlePt(const edm4hep::ReconstructedParticle& p1,
                               const edm4hep::ReconstructedParticle& p2) const {
  // transverse momentum components
  double px_tot = p1.getMomentum().x + p2.getMomentum().x;
  double py_tot = p1.getMomentum().y + p2.getMomentum().y;

  return std::sqrt(px_tot * px_tot + py_tot * py_tot);
}

// set up root tree branches
void SelectEvents::setupBranches() {
  tree->Branch("runNumber", &runNumber, "runNumber/I");
  tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");
  tree->Branch("weight", &weight, "weight/F");
  tree->Branch("sqs_eff", &sqs_eff, "sqs_eff/F");
  tree->Branch("channelType", &channelType, "channelType/I");
  tree->Branch("processID", &proc_id, "processID/I");
  tree->Branch("passElectronPreSelection", &passElectronPreSelection, "passElectronPreSelection/B");
  tree->Branch("passMuonPreSelection", &passMuonPreSelection, "passMuonPreSelection/B");
  tree->Branch("passJetPreSelection", &passJetPreSelection, "passJetPreSelection/B");

  // Jets
  tree->Branch("nJets", &nJets, "nJets/I");
  tree->Branch("jet_pt", jet_pt, "jet_pt[nJets]/F");
  tree->Branch("jet_eta", jet_eta, "jet_eta[nJets]/F");
  tree->Branch("jet_phi", jet_phi, "jet_phi[nJets]/F");
  tree->Branch("jet_e", jet_e, "jet_e[nJets]/F");
  tree->Branch("jet_m", jet_m, "jet_m[nJets]/F");

  // electrons
  tree->Branch("nElectrons", &nElectrons, "nElectrons/I");
  tree->Branch("ele_pt", ele_pt, "ele_pt[nElectrons]/F");
  tree->Branch("ele_eta", ele_eta, "ele_eta[nElectrons]/F");
  tree->Branch("ele_phi", ele_phi, "ele_phi[nElectrons]/F");
  tree->Branch("ele_e", ele_e, "ele_e[nElectrons]/F");
  tree->Branch("ele_charge", ele_charge, "ele_charge[nElectrons]/I");
  tree->Branch("di_electron_mass", &di_electron_mass, "di_electron_mass/F");
  tree->Branch("di_electron_pt", &di_electron_pt, "di_electron_pt/F");

  // muons
  tree->Branch("nMuons", &nMuons, "nMuons/I");
  tree->Branch("muon_pt", muon_pt, "muon_pt[nMuons]/F");
  tree->Branch("muon_eta", muon_eta, "muon_eta[nMuons]/F");
  tree->Branch("muon_phi", muon_phi, "muon_phi[nMuons]/F");
  tree->Branch("muon_e", muon_e, "muon_e[nMuons]/F");
  tree->Branch("muon_charge", muon_charge, "muon_charge[nMuons]/I");
  tree->Branch("di_muon_mass", &di_muon_mass, "di_muon_mass/F");
  tree->Branch("di_muon_pt", &di_muon_pt, "di_muon_pt/F");

  tree->Branch("MET", &MET, "MET/F");
  tree->Branch("MET_phi", &MET_phi, "MET_phi/F");
  tree->Branch("recoil_mass", &recoil_mass, "recoil_mass/F");
  tree->Branch("recoil_pt", &recoil_pt, "recoil_pt/F");
  tree->Branch("recoil_cosTheta", &recoil_cosTheta, "recoil_cosTheta/F");
  // kinematics
  tree->Branch("acoplanarity", &acoplanarity, "acoplanarity/F");
  tree->Branch("acollinearity", &acollinearity, "acollinearity/F");
  tree->Branch("deltaPhi", &deltaPhi, "deltaPhi/F");
  tree->Branch("deltaR", &deltaR, "deltaR/F");
  tree->Branch("visible_mass", &visible_mass, "visible_mass/F");
  tree->Branch("visible_energy", &visible_energy, "visible_energy/F");
  tree->Branch("visible_pt", &visible_pt, "visible_pt/F");

  // MC truth
  tree->Branch("nMCParticles", &nMCParticles, "nMCParticles/I");
  tree->Branch("mc_pt", &mc_pt);
  tree->Branch("mc_eta", &mc_eta);
  tree->Branch("mc_phi", &mc_phi);
  tree->Branch("mc_e", &mc_e);
  tree->Branch("mc_pdgId", &mc_pdgId);
  tree->Branch("mc_status", &mc_status);
  tree->Branch("mc_motherPdgId", &mc_motherPdgId);
}

void SelectEvents::addMCParticle(const TLorentzVector& p4, int pdgId, int status, int motherPdgId) const {
  addMCParticle(p4.Pt(), p4.Eta(), p4.Phi(), p4.E(), pdgId, status, motherPdgId);
}

void SelectEvents::addMCParticle(float pt, float eta, float phi, float e, int pdgId, int status, int motherPdgId) const {
  mc_pt.push_back(pt);
  mc_eta.push_back(eta);
  mc_phi.push_back(phi);
  mc_e.push_back(e);
  mc_pdgId.push_back(pdgId);
  mc_status.push_back(status);
  mc_motherPdgId.push_back(motherPdgId);
  nMCParticles = mc_pt.size();
}

void SelectEvents::fillEvent() const {
  tree->Fill();
  // Clear MC truth vectors
  mc_pt.clear();
  mc_eta.clear();
  mc_phi.clear();
  mc_e.clear();
  mc_pdgId.clear();
  mc_status.clear();
  mc_motherPdgId.clear();
}
