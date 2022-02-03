#include "Geometry/HGCalCommonData/interface/HGCalCellUV.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <vector>

HGCalCellUV::HGCalCellUV(double waferSize, double separation, int32_t nFine, int32_t nCoarse) {
  HGCalCell hgcalcell(waferSize, nFine, nCoarse);
  ncell_[0] = nFine;
  ncell_[1] = nCoarse;
  for (int k = 0; k < 2; ++k) {
    cellX_[k] = waferSize / (3 * ncell_[k]);
    cellY_[k] = 0.5 * sqrt3_ * cellX_[k];
    cellXTotal_[k] = (waferSize + separation) / (3 * ncell_[k]);
    cellY_[k] = 0.5 * sqrt3_ * cellXTotal_[k];
  }

  // Fill up the placement vectors
  for (int placement = 0; placement < HGCalCell::cellPlacementTotal; ++placement) {
    // Fine cells
    for (int iu = 0; iu < 2 * ncell_[0]; ++iu) {
      for (int iv = 0; iv < 2 * ncell_[0]; ++iv) {
	int u = (placement < HGCalCell::cellPlacementExtra) ? iv : iu;
	int v = (placement < HGCalCell::cellPlacementExtra) ? iu : iv;
	if (((v - u) < ncell_[0]) && (u - v) <= ncell_[0]) {
	  cellPosFine_[placement][std::pair<int, int>(u, v)] = hgcalcell.HGCalCellUV2XY1(u, v, placement, 0);
	}
      }
    }
    // Coarse cells
    for (int iu = 0; iu < 2 * ncell_[1]; ++iu) {
      for (int iv = 0; iv < 2 * ncell_[1]; ++iv) {
	int u = (placement < HGCalCell::cellPlacementExtra) ? iv : iu;
	int v = (placement < HGCalCell::cellPlacementExtra) ? iu : iv;
	if (((v - u) < ncell_[1]) && (u - v) <= ncell_[1]) {
	  cellPosCoarse_[placement][std::pair<int, int>(u, v)] = hgcalcell.HGCalCellUV2XY1(u, v, placement, 1);
	}
      }
    }
  }
}


std::pair<int32_t, int32_t> HGCalCellUV::HGCalCellUVFromXY1(double xloc, double yloc, int32_t placement, int32_t type, bool extend, bool debug) {
  double y = yloc;
  double x = (placement>=6) ? xloc: -1*xloc;
  const std::vector<double> fcos = {0.5, 1.0, 0.5, -0.5, -1.0, -0.5};
  const std::vector<double> fsin = {-sqrt3By2_, 0.0, sqrt3By2_, sqrt3By2_, 0.0, -sqrt3By2_};
  x = x * fcos[placement%6] + y * fsin[placement%6];
  y = x * fsin[placement%6] + y * fcos[placement%6];
  
  int32_t u(-100), v(-100);
  double r2 = 0.5 * sqrt3_;
  double R2 = 2 * r2 / sqrt3_;
  int N_ = (type!=0) ? ncell_[1]: ncell_[0]; 
  double c1 = y + (x / sqrt3_);
  double c2 = y - (x / sqrt3_);
  if ((x > r2) || (x < -1*r2) || (c1 > R2) || (c1 < -1*R2) || (c2 > R2) || (c2 < -1*R2)){
    return std::make_pair(u, v);
  }
  else{
    double R = waferSize / (3 * N_);
    double r = sqrt3By2_ * R;
    int l1 = std::floor((y/r) + N_ -1);
    int l2 = std::floor((0.5*y + 0.5*x/sqrt3_)/r + N_ -4/3);
    int l3 = std::floor((x/sqrt3_)/r + N_ -4/3);
    int u1 = std::ceil((y/r) + N_ +1);
    int u2 = std::ceil((0.5*y + 0.5*x/sqrt3_)/r + N_ -2/3);
    int u3 = std::ceil((x/sqrt3_)/r + N_);
    //std::cout<< l1 << "   " << l2 << "   " << l3 << "   " << u1 << "   " << u2 << "   " << u3 << std::endl;
    for(int ui = l2 +1; ui<u2; ui++){
      for(int vi = l3 +1; vi<u3; vi++){
        int c3 = 2*ui - vi;
        //std::cout<< c3 << std::endl;
        if((c3<u1)||(c3>l1)){
          u = ui;
          v = vi;   
        }
      }
    }
  return std::make_pair(u, v);
  }

}



std::pair<int, int> HGCalCellUV::HGCalCellUVFromXY2(double xloc, double yloc, int32_t placement, int32_t type, bool extend, bool debug) {
  if (type != 0)
    return HGCalCellUVFromXY2(xloc, yloc, ncell_[1], cellX_[1], cellY_[1], cellXTotal_[1], cellY_[1], cellPosCoarse_[placement], extend, debug);
  else
    return HGCalCellUVFromXY2(xloc, yloc, ncell_[0], cellX_[0], cellY_[0], cellXTotal_[0], cellY_[0], cellPosFine_[placement], extend, debug);
}

std::pair<int, int> HGCalCellUV::HGCalCellUVFromXY2(double xloc, double yloc, int n, double cellX, double cellY, double cellXTotal, double cellYTotal, std::map<std::pair<int, int>, std::pair<double, double> >& cellPos, bool extend, bool debug) {
  std::pair<int, int> uv = std::make_pair(-1, -1);
  std::map<std::pair<int, int>, std::pair<double, double> >::const_iterator itr;
  for (itr = cellPos.begin(); itr != cellPos.end(); ++itr) {
    double delX = std::abs(xloc - (itr->second).first);
    double delY = std::abs(yloc - (itr->second).second);
    if ((delX < cellX) && (delY < cellY)) {
      if ((delX < (0.5 * cellX)) || (delY < (2.0 * cellY - sqrt3_ * delX))) {
	uv = itr->first;
	break;
      }
    }
  }
  if ((uv.first < 0) && extend) {
    for (itr = cellPos.begin(); itr != cellPos.end(); ++itr) {
      double delX = std::abs(xloc - (itr->second).first);
      double delY = std::abs(yloc - (itr->second).second);
      if ((delX < cellXTotal) && (delY < cellYTotal)) {
	if ((delX < (0.5 * cellXTotal)) || (delY < (2.0 * cellYTotal - sqrt3_ * delX))) {
	  uv = itr->first;
	  break;
	}
      }
    }
  }
  if (debug) edm::LogVerbatim("HGCalGeom") << "HGCalCellUVFromXY2: Input " << xloc << ":" << yloc << ":" << extend << " Output " << uv.first << ":" << uv.second;
  return uv;
}
