#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using LorentzVectorXYZ = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;
using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecS = ROOT::VecOps::RVec<size_t>;
using RVecUC = ROOT::VecOps::RVec<UChar_t>;
using RVecF = ROOT::VecOps::RVec<float>;
using RVecB = ROOT::VecOps::RVec<bool>;
using RVecVecI = ROOT::VecOps::RVec<RVecI>;
using RVecLV = ROOT::VecOps::RVec<LorentzVectorM>;
using RVecSetInt = ROOT::VecOps::RVec<std::set<int>>;

ROOT::VecOps::RVec<LorentzVectorM> yumengGetP4(
    const ROOT::VecOps::RVec<double>& pt,
    const ROOT::VecOps::RVec<double>& eta,
    const ROOT::VecOps::RVec<double>& phi,
    const ROOT::VecOps::RVec<double>& mass) {

    ROOT::VecOps::RVec<LorentzVectorM> p4(pt.size());
    
    for (size_t i = 0; i < pt.size(); i++) {
        p4[i] = LorentzVectorM(
            static_cast<double>(pt[i]),
            static_cast<double>(eta[i]),
            static_cast<double>(phi[i]),
            static_cast<double>(mass[i])
        );
    }   
    return p4;
}

float FindMatchingdR(const LorentzVectorM& target_p4, const RVecLV& ref_p4,const float deltaR_thr){
  double deltaR_min = deltaR_thr;
  int current_idx = -1;
  for(int refIdx =0; refIdx<ref_p4.size(); refIdx++){
    auto dR_targetRef= ROOT::Math::VectorUtil::DeltaR(target_p4, ref_p4.at(refIdx));
    if ( dR_targetRef < deltaR_min ) {
      deltaR_min = dR_targetRef ;
      current_idx = refIdx;
    }
  }
  return deltaR_min;
}

std::set<size_t> FindMatchingSet(const LorentzVectorM& target_p4, const RVecLV& ref_p4, const float deltaR_thr){
  std::set<size_t> matched_idx;
  for(size_t refIdx =0; refIdx<ref_p4.size(); ++refIdx){
    auto dR_targetRef= ROOT::Math::VectorUtil::DeltaR(target_p4, ref_p4.at(refIdx));
    if ( dR_targetRef < deltaR_thr ) {
      matched_idx.insert(refIdx);
    }
  }
  return matched_idx;
}
