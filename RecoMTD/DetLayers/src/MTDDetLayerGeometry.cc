/** \file
 *
 *  \author L. Gray - FNAL
 *
 */

#include <RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h>
#include <RecoMTD/DetLayers/interface/ETLDetLayerGeometryBuilder.h>
#include <RecoMTD/DetLayers/interface/BTLDetLayerGeometryBuilder.h>

#include <FWCore/Utilities/interface/Exception.h>
#include <TrackingTools/DetLayers/interface/DetLayer.h>
#include <DataFormats/ForwardDetId/interface/BTLDetId.h>
#include <DataFormats/ForwardDetId/interface/ETLDetId.h>

#include <Utilities/General/interface/precomputed_value_sort.h>
#include <DataFormats/GeometrySurface/interface/GeometricSorting.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <algorithm>

using namespace std;
using namespace geomsort;
using namespace edm;

MTDDetLayerGeometry::MTDDetLayerGeometry() {}

MTDDetLayerGeometry::~MTDDetLayerGeometry() {}

void MTDDetLayerGeometry::buildLayers(const MTDGeometry* geo, const MTDTopology* mtopo) {
  bool abort(false);
  if (geo == nullptr) {
    LogError("MTDDetLayers") << "No MTD geometry is available.";
    abort = true;
  }
  if (mtopo == nullptr) {
    LogError("MTDDetLayers") << "No MTD topology  is available.";
    abort = true;
  }
  if (abort) {
    throw cms::Exception("MTDDetLayers") << "No complete MTD geometry available, aborting.";
  }

  // Build BTL layers
  this->addBTLLayers(BTLDetLayerGeometryBuilder::buildLayers(*geo, *mtopo));
  // Build ETL layers, depends on the scenario
  this->addETLLayers(ETLDetLayerGeometryBuilder::buildLayers(*geo, *mtopo));
}

void MTDDetLayerGeometry::addETLLayers(const pair<vector<DetLayer*>, vector<DetLayer*> >& etllayers) {
  for (auto const it : etllayers.first) {
    etlLayers_fw.push_back(it);
    allForward.push_back(it);

    detLayersMap[makeDetLayerId(it)] = it;
  }

  for (auto const it : etllayers.second) {
    etlLayers_bk.push_back(it);
    allBackward.push_back(it);

    detLayersMap[makeDetLayerId(it)] = it;
  }
}

void MTDDetLayerGeometry::addBTLLayers(const vector<DetLayer*>& dtlayers) {
  for (auto const it : dtlayers) {
    btlLayers.push_back(it);
    allBarrel.push_back(it);

    detLayersMap[makeDetLayerId(it)] = it;
  }
}

DetId MTDDetLayerGeometry::makeDetLayerId(const DetLayer* detLayer) const {
  if (detLayer->subDetector() == GeomDetEnumerators::TimingEndcap) {
    ETLDetId id(detLayer->basicComponents().front()->geographicalId().rawId());
    return ETLDetId(id.mtdSide(), 0, 0, 0, 0);  // Constructor of new geometry is compatible with prev8
  } else if (detLayer->subDetector() == GeomDetEnumerators::TimingBarrel) {
    BTLDetId id(detLayer->basicComponents().front()->geographicalId().rawId());
    return BTLDetId(id.mtdSide(), 0, 0, 0, 0, 0);
  } else
    throw cms::Exception("InvalidModuleIdentification");  // << detLayer->module();
}

const vector<const DetLayer*>& MTDDetLayerGeometry::allBarrelLayers() const { return allBarrel; }

const vector<const DetLayer*>& MTDDetLayerGeometry::allEndcapLayers() const { return allEndcap; }

const vector<const DetLayer*>& MTDDetLayerGeometry::allForwardLayers() const { return allForward; }

const vector<const DetLayer*>& MTDDetLayerGeometry::allBackwardLayers() const { return allBackward; }

const vector<const DetLayer*>& MTDDetLayerGeometry::allBTLLayers() const { return btlLayers; }

const vector<const DetLayer*>& MTDDetLayerGeometry::allETLLayers() const { return etlLayers_all; }

const vector<const DetLayer*>& MTDDetLayerGeometry::allLayers() const { return allDetLayers; }

////////////////////////////////////////////////////

const DetLayer* MTDDetLayerGeometry::idToLayer(const DetId& id) const {
  DetId idout;
  MTDDetId detId(id);

  if (detId.mtdSubDetector() == 2) {  // 2 is ETL
    ETLDetId etlId(detId.rawId());
    idout = ETLDetId(etlId.mtdSide(), 0, 0, 0, 0);
  } else if (detId.mtdSubDetector() == 1) {  // 1 is BTL
    BTLDetId btlId(detId.rawId());
    idout = BTLDetId(btlId.mtdSide(), 0, 0, 0, 0, 0);
  } else
    throw cms::Exception("InvalidSubdetId") << detId.subdetId();

  std::map<DetId, const DetLayer*>::const_iterator layer = detLayersMap.find(idout);
  if (layer == detLayersMap.end())
    return nullptr;
  return layer->second;
}

// Quick way to sort barrel det layers by increasing R,
// do not abuse!
#include <TrackingTools/DetLayers/interface/BarrelDetLayer.h>
struct ExtractBarrelDetLayerR {
  typedef Surface::Scalar result_type;
  result_type operator()(const DetLayer* p) const {
    const BarrelDetLayer* bdl = dynamic_cast<const BarrelDetLayer*>(p);
    if (bdl)
      return bdl->specificSurface().radius();
    else
      return -1.;
  }
};

void MTDDetLayerGeometry::sortLayers() {
  // The following are filled inside-out, no need to re-sort
  // precomputed_value_sort(dtLayers.begin(), dtLayers.end(),ExtractR<DetLayer,float>());
  // precomputed_value_sort(cscLayers_fw.begin(), cscLayers_fw.end(),ExtractAbsZ<DetLayer,float>());
  // precomputed_value_sort(cscLayers_bk.begin(), cscLayers_bk.end(),ExtractAbsZ<DetLayer,float>());
  // precomputed_value_sort(rpcLayers_fw.begin(), rpcLayers_fw.end(),ExtractAbsZ<DetLayer,float>());
  // precomputed_value_sort(rpcLayers_bk.begin(), rpcLayers_bk.end(),ExtractAbsZ<DetLayer,float>());
  // precomputed_value_sort(rpcLayers_barrel.begin(), rpcLayers_barrel.end(), ExtractR<DetLayer,float>());

  // Sort these inside-out
  precomputed_value_sort(allBarrel.begin(), allBarrel.end(), ExtractBarrelDetLayerR());
  precomputed_value_sort(allBackward.begin(), allBackward.end(), ExtractAbsZ<DetLayer, float>());
  precomputed_value_sort(allForward.begin(), allForward.end(), ExtractAbsZ<DetLayer, float>());

  // Build more complicated vectors with correct sorting

  //etlLayers_all: from -Z to +Z
  etlLayers_all.reserve(etlLayers_bk.size() + etlLayers_fw.size());
  std::copy(etlLayers_bk.begin(), etlLayers_bk.end(), back_inserter(etlLayers_all));
  std::reverse(etlLayers_all.begin(), etlLayers_all.end());
  std::copy(etlLayers_fw.begin(), etlLayers_fw.end(), back_inserter(etlLayers_all));

  // allEndcap: order is  all bw, all fw
  allEndcap.reserve(allBackward.size() + allForward.size());
  std::copy(allBackward.begin(), allBackward.end(), back_inserter(allEndcap));
  std::reverse(allEndcap.begin(), allEndcap.end());
  std::copy(allForward.begin(), allForward.end(), back_inserter(allEndcap));

  // allDetLayers: order is  all bw, all barrel, all fw
  allDetLayers.reserve(allBackward.size() + allBarrel.size() + allForward.size());
  std::copy(allBackward.begin(), allBackward.end(), back_inserter(allDetLayers));
  std::reverse(allDetLayers.begin(), allDetLayers.end());
  std::copy(allBarrel.begin(), allBarrel.end(), back_inserter(allDetLayers));
  std::copy(allForward.begin(), allForward.end(), back_inserter(allDetLayers));

  // number layers
  int sq = 0;
  for (auto l : allDetLayers)
    (*const_cast<DetLayer*>(l)).setSeqNum(sq++);
}
