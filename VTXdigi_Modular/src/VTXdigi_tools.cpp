// VTXdigi_Modular/src/VTXdigi_tools.cpp
#include "VTXdigi_tools.h"
#include <DD4hep/Objects.h>
#include <DD4hep/VolumeManager.h>
#include <edm4hep/Vector3d.h>

namespace VTXdigi_tools {

SimHitWrapper::SimHitWrapper(
  edm4hep::SimTrackerHit simTrackerHit, dd4hep::DDSegmentation::VolumeID volumeID,
  const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder,
  const dd4hep::VolumeManager& volumeManager,
  const std::unique_ptr<dd4hep::rec::CellIDPositionConverter>& cellIDPositionConverter)
    : m_simTrackerHit(simTrackerHit), m_volumeID(volumeID) {

  m_charge = static_cast<float>(m_simTrackerHit.getEDep() * (dd4hep::GeV / dd4hep::keV) * kChargePerkeV); // convert energy deposit (in keV) to number of electrons
  m_layerNumber = GetLayer(m_volumeID, cellIdDecoder);

  // check if the simHit was caused by a primary, secondary or delta particle
  if ( m_simTrackerHit.isProducedBySecondary() ) {
    // if ddsim dropped the MCParticle that caused this simHit, we assume it was a delta ray
    // ddsim drops MCParticles below a certain energy cut to save computing cost and disk space.
    m_mcParticleLevel = MCParticleLevel::Delta;
  }
  else {
    // check if the MCPArticle was created by the generator
    const int32_t simulatorStatus = m_simTrackerHit.getParticle().getSimulatorStatus();
    const int32_t mask = 1 << edm4hep::MCParticle::BITCreatedInSimulation; // should be bit 30
    const bool causedByPrimary = (simulatorStatus & mask) == 0; // bit is not set -> created in generator
    if ( causedByPrimary )
      m_mcParticleLevel = MCParticleLevel::Primary;
    else
    {
      // now check if the MCParticle prod. vertex lies outside this sensors volume (by comparing volumeIDs)
      const edm4hep::Vector3d prodVertex_temp = m_simTrackerHit.getParticle().getVertex();
      const dd4hep::Position prodVertex = 0.1 * dd4hep::Position(prodVertex_temp.x, prodVertex_temp.y, prodVertex_temp.z); // convert edm4hep's mm -> dd4hep's cm
      const dd4hep::DDSegmentation::CellID prodVertex_cellID = cellIDPositionConverter->cellID(prodVertex); // returns 0 if the position is outside of any sensitive volume

      // convert cellID to volumeID - see comment in VTXdigi_Modular::GetVolumeID())
      dd4hep::DDSegmentation::CellID prodVertex_volumeID;
      if (prodVertex_cellID == 0)
        prodVertex_volumeID = 0; // lookupContext(cellID=0) crashes (because cellID 0 does not exist)
      else
        prodVertex_volumeID = volumeManager.lookupContext(prodVertex_cellID)->element.volumeID();

      if (prodVertex_volumeID != m_volumeID) {
        // the MCParticle was created outside of this sensor's sensitive volume
        m_mcParticleLevel = MCParticleLevel::Secondary;
      }
      else {
        m_mcParticleLevel = MCParticleLevel::Delta;
      }
    }
  }
}

void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept {
  std::swap(a.m_simTrackerHit, b.m_simTrackerHit);
  std::swap(a.m_volumeID, b.m_volumeID);
  std::swap(a.m_charge, b.m_charge);
  std::swap(a.m_layerNumber, b.m_layerNumber);
  std::swap(a.m_truthPos, b.m_truthPos);
  std::swap(a.m_mcParticleLevel, b.m_mcParticleLevel);
} // swap(Hit&, Hit&)

// SimulatorStatus bits (see https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_m_c_particle.html)
// 29 : "Backscatter",
// 30 : "CreatedInSimulation",
// 26 : "DecayedInCalorimeter",
// 27 : "DecayedInTracker",
// 22 : "HandledInFastSim",
// 25 : "LeftWorld",
// 23 : "Overlay",
// 24 : "Stopped",
// 28 : "VertexIsNotEndpointOfParent",

/* -- helpers -- */

dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec) {
  return dd4hep::rec::Vector3D(vec.x, vec.y, vec.z);
}
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3f vec) {
  return dd4hep::rec::Vector3D(static_cast<double>(vec.x), static_cast<double>(vec.y), static_cast<double>(vec.z));
}
edm4hep::Vector3d ConvertVector(dd4hep::rec::Vector3D vec) {
  return edm4hep::Vector3d(vec.x(), vec.y(), vec.z());
}

TGeoHMatrix ComputeSensorTrafoMatrix(const dd4hep::DDSegmentation::VolumeID& volumeID, const dd4hep::VolumeManager& volumeManager, const TGeoRotation& sensorNormalRotation) {
  TGeoHMatrix M = volumeManager.lookupDetElement(volumeID).nominal().worldTransformation();

  /* rotate the local coordinate system st. sensor U is (1,0,0), V is (0,1,0) and normal vector is (0,0,1) */
  M.Multiply(sensorNormalRotation);

  /* rotation is unitless, but need to convert translation from cm to mm (dd4hep::mm = 0.1) */
  double* transl = M.GetTranslation();
  transl[0] = transl[0] / dd4hep::mm;
  transl[1] = transl[1] / dd4hep::mm;
  transl[2] = transl[2] / dd4hep::mm;
  M.SetTranslation(transl);

  return M;
}

dd4hep::rec::Vector3D Trafo_global_local(const dd4hep::rec::Vector3D& global, const TGeoHMatrix& M) {
  double local[3];
  M.MasterToLocal(global, local);
  return dd4hep::rec::Vector3D(local[0], local[1], local[2]);
}

dd4hep::rec::Vector3D Trafo_local_global(const dd4hep::rec::Vector3D& local, const TGeoHMatrix& M) {
  double global[3];
  M.LocalToMaster(local, global);
  return dd4hep::rec::Vector3D(global[0], global[1], global[2]);
}

std::array<float, 2> Trafo_local_pixIndexCoords(const dd4hep::rec::Vector3D& local, const std::array<float, 2> pixelPitch, const std::array<size_t, 2> pixelCount) {
  const std::array<float, 2> local_2d = {static_cast<float>(local.x()), static_cast<float>(local.y())};
  std::array<float, 2> pixIndex;
  for (size_t axis = 0; axis < 2; ++axis) {
    const float halfLength = 0.5 * pixelPitch[axis] * pixelCount[axis];
    if (local_2d[axis] < -halfLength || local_2d[axis] > halfLength) {
      throw std::runtime_error("VTXdigi_tools::ComputePixelIndexCoords(): position is out of sensor bounds");
    }
    pixIndex[axis] = (local_2d[axis] + halfLength) / pixelPitch[axis] - 0.5f; // shift from [-halfLength, halfLength] to [-0.5, pixelCount - 0.5]
  }
  return pixIndex;
}

dd4hep::rec::Vector3D Trafo_pixIndexCoords_local(const std::array<float, 2>& pixIndexCoords, const float w, const std::array<float, 2> pixelPitch, const std::array<size_t, 2> pixelCount) {
  std::array<float, 2> local;
  for (size_t axis = 0; axis < 2; ++axis) {
    const float halfLength = 0.5 * pixelPitch[axis] * pixelCount[axis];
    local[axis] = (pixIndexCoords[axis] + 0.5f) * pixelPitch[axis] - halfLength; // shift from [-0.5, pixelCount - 0.5] to [-halfLength, halfLength]
  }
  return dd4hep::rec::Vector3D(local[0], local[1], w);
}

int GetLayer(const dd4hep::DDSegmentation::VolumeID& volumeID, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) {
  return static_cast<int>(cellIdDecoder->get(volumeID, "layer"));
}



/* -- Binning things -- */

int ComputeBinIndex(float x, float binX0, float binWidth, int binN) {
  /** Get the bin index for a given x value
   *  binX0 is the lower edge of the first bin
   *  binWidth is the width of the bins
   *  binN is the number of bins
   *  return -1 if x is out of range
   */

  if (binN <= 0) throw std::runtime_error("VTXdigi_tools::ComputeBinIndex(): binN must be positive");
  if (binWidth <= 0.0) throw std::runtime_error("VTXdigi_tools::ComputeBinIndex(): binWidth must be positive");

  float relativePos = (x - binX0) / binWidth; // shift to [0, binN]
  if (relativePos < 0.0f || relativePos > static_cast<float>(binN))
    return -1;
  if (relativePos == static_cast<float>(binN))
    return binN - 1; // include upper edge in last bin (makes sense for pixels)
  return static_cast<int>(relativePos);
} // ComputeBinIndex()

float ComputeBinCenter(int i, float binX0, float binWidth) {
  if (i < 0) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): bin index must be non-negative");
  if (binWidth <= 0.0) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): binWidth must be positive");

  return binX0 + (static_cast<float>(i) + 0.5f) * binWidth; // add 0.5*binWidth to shift from lower edge to center
} // ComputeBinCenter()
float ComputeBinCenter(int i, float binX0, float binX1, int binN) {
  if (binN <= 0) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): binN must be positive");
  if (binX1 <= binX0) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): binX1 must be greater than binX0");
  if (i < 0 || i >= binN) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): bin index out of bounds");

  const float binWidth = (binX1 - binX0) / static_cast<float>(binN);
  return ComputeBinCenter(i, binX0, binWidth);
} // ComputeBinCenter()

std::array<int, 2> Trafo_local_pixIndex(const dd4hep::rec::Vector3D& pos, const std::array<float, 2> pixelPitch, const std::array<size_t, 2> pixelCount) {
  const float length_u_half = 0.5 * pixelPitch[0] * pixelCount[0];
  int i_u = ComputeBinIndex(
    pos.x(),
    -length_u_half,
    pixelPitch[0],
    pixelCount[0]);

  const float length_v_half = 0.5 * pixelPitch[1] * pixelCount[1];
  int i_v = ComputeBinIndex(
    pos.y(),
    -length_v_half,
    pixelPitch[1],
    pixelCount[1]);

  return {i_u, i_v};
} // Trafo_local_pixIndex()


std::array<int, 3> Trafo_local_inpixIndex(const dd4hep::rec::Vector3D& pos, const std::array<int, 3>& binCount, const std::array<float, 2>& pixelPitch, const std::array<float, 3>& activeVolumeDimensions) {
  std::array<int, 3> indices;

  const float posShifted_u = pos.x() + 0.5 * activeVolumeDimensions[0]; // shift to [0, length_u]
  if (posShifted_u < 0.0 || posShifted_u > activeVolumeDimensions[0]) {
    indices[0] = -1; // out of bounds
  }
  else {
    float posInPixel_u = std::fmod(posShifted_u,  pixelPitch[0]);
    if (posInPixel_u < 0.0) posInPixel_u +=  pixelPitch[0]; // ensure positive remainder
    indices[0] = ComputeBinIndex(posInPixel_u, 0.0,  pixelPitch[0] / binCount[0], binCount[0]);
  }

  const float posShifted_v = pos.y() + 0.5 * activeVolumeDimensions[1];
  if (posShifted_v < 0.0 || posShifted_v > activeVolumeDimensions[1]) {
    indices[1] = -1; // out of bounds
  }
  else {
    float posInPixel_v = std::fmod(posShifted_v, pixelPitch[1]);
    if (posInPixel_v < 0.0) posInPixel_v += pixelPitch[1];
    indices[1] = ComputeBinIndex(posInPixel_v, 0.0, pixelPitch[1] / binCount[1], binCount[1]);
  }

  // vertical (w) binning: shift to [0, thickness]
  const float posShifted_w = pos.z() + 0.5 * activeVolumeDimensions[2];
  indices[2] = ComputeBinIndex(posShifted_w, 0.0, activeVolumeDimensions[2] / binCount[2], binCount[2]); // no fmod, so out-of-bounds is caught

  return indices;
} // Trafo_local_inpixIndex()

dd4hep::rec::Vector3D Trafo_pixIndex_local(const std::array<int, 2> pixelIndex, const std::array<float, 2> sensorLength,  const std::array<float, 2> pixelPitch, float depletedRegionDepthCenter) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */

  float u = (static_cast<float>(pixelIndex[0]) + 0.5f) * pixelPitch[0] - 0.5f * sensorLength[0]; // in mm
  float v = (static_cast<float>(pixelIndex[1]) + 0.5f) * pixelPitch[1] - 0.5f * sensorLength[1];
  float w = depletedRegionDepthCenter;

  return dd4hep::rec::Vector3D(u, v, w);
}

dd4hep::rec::Vector3D Trafo_pixIndex_local(const std::array<int, 2> pixelIndex, const std::array<float, 2> sensorLength, const std::array<float, 2> pixelPitch) {
  return Trafo_pixIndex_local(pixelIndex, sensorLength, pixelPitch, 0.f);
}

dd4hep::rec::Vector3D Trafo_pixIndex_local(const std::array<float, 2> index, const std::array<float, 2> sensorLength,  const std::array<float, 2> pixelPitch, float depletedRegionDepthCenter) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */

  float u = (index[0] + 0.5f) * pixelPitch[0] - 0.5f * sensorLength[0]; // in mm. Add 0.5*pixelPitch to shift from pixel edge to center, since index 0 is defined as the center of the pixel.
  float v = (index[1] + 0.5f) * pixelPitch[1] - 0.5f * sensorLength[1];
  float w = depletedRegionDepthCenter;

  return dd4hep::rec::Vector3D(u, v, w);
}

dd4hep::rec::Vector3D Trafo_pixIndex_local(const std::array<float, 2> index, const std::array<float, 2> sensorLength,  const std::array<float, 2> pixelPitch) {
  return Trafo_pixIndex_local(index, sensorLength, pixelPitch, 0.f);
}

/* -- HitMap -- */

HitMap::HitMap(std::array<size_t, 2> pixelCount) : m_pixCount(pixelCount) {
  const int inverseOccupancy = 2000; // assume occupancy, 5e-4 is quite conservative for Z-run
  m_pixels.reserve(pixelCount[0] * pixelCount[1] / inverseOccupancy); // avoid too many reallocations
}

void HitMap::FillCharge(std::array<int, 2> i_uv, float charge, const SimHitWrapper& simHitWrapper) {
  if (charge < 1.e-6f)
    return; // skip very small charge additions for performance (this is NECESSARY to skip in-pix bins with weight ~0)
  if (_OutOfBounds(i_uv)) [[unlikely]]
    throw std::runtime_error("HitMap::FillCharge: pixel i_u or i_v ( " + std::to_string(i_uv[0]) + ", " + std::to_string(i_uv[1]) + ") out of range");

  auto [iter, inserted] = m_pixels.try_emplace(i_uv, Pixel(i_uv));
  iter->second.charge += charge;
  iter->second.simHits.insert(&simHitWrapper);
}

void HitMap::ApplyChargeSmearing(const Rndm::Numbers& rndm_charge) {
  auto hitIter = m_pixels.begin();
  while (hitIter != m_pixels.end()) {
    hitIter->second.charge = std::max(hitIter->second.charge + static_cast<float>(rndm_charge()), 0.f); // don't allow negative charge after smearing
    ++hitIter;
  }
}

void HitMap::ApplyThreshold(const float threshold, const Rndm::Numbers* rndm_threshold) {
  auto hitIter = m_pixels.begin();
  while (hitIter != m_pixels.end()) {
    // optionally disperse the threshold per pixel (drawn per event per sensor per pixel)
    const float pixThreshold = rndm_threshold ? threshold + static_cast<float>((*rndm_threshold)()) : threshold;
    if (hitIter->second.charge < pixThreshold)
      hitIter = m_pixels.erase(hitIter); // erase returns the iterator to the next element, so this is safe to do while iterating
    else
      ++hitIter;
  }
}

float HitMap::GetCharge(std::array<int, 2> i_uv) const {
  if (_OutOfBounds(i_uv)) [[unlikely]] {
    throw std::runtime_error("HitMap::GetCharge: pixel i_u or i_v ( " + std::to_string(i_uv[0]) + ", " + std::to_string(i_uv[1]) + ") out of range");
  }
  auto it = m_pixels.find(i_uv);
  if (it == m_pixels.end())
    return 0.f; // if pixel not found, charge is 0
  return it->second.charge;
}

float HitMap::GetTotalCharge() const {
  float totalCharge = 0.f;
  for (const auto& [i_uv, pixHit] : m_pixels) {
    totalCharge += pixHit.charge;
  }
  return totalCharge;
}

inline bool HitMap::_OutOfBounds(std::array<int, 2> i_uv) const {
  return (
    i_uv[0] < 0
    || i_uv[0] >= static_cast<int>(m_pixCount[0])
    || i_uv[1] < 0
    || i_uv[1] >= static_cast<int>(m_pixCount[1])
  );
}

/* -- Eta function -- */

EtaDistribution::EtaDistribution(std::array<std::vector<std::pair<float, float>>, 2> points) : m_points(points){
  for (int i_axis=0; i_axis < 2; ++i_axis) {
    if (m_points[i_axis].empty())
      throw std::runtime_error("EtaDistribution: no distribution points provided for axis " + std::to_string(i_axis));
    if (m_points[i_axis].size() < 2)
      throw std::runtime_error("EtaDistribution: need at least 2 distribution points for axis " + std::to_string(i_axis));

    // points need to be delivered sorted by biased position (first element of pair) for linear interpolation to work.
    for (size_t i_point=0; i_point < m_points[i_axis].size()-1; ++i_point) {
      if (m_points[i_axis][i_point].first <0 || m_points[i_axis][i_point].first > 1.0)
        throw std::runtime_error("EtaDistribution: biased position for axis " + std::to_string(i_axis) + " is out of bounds [0,1] for point " + std::to_string(i_point) + " (biased pos. " + std::to_string(m_points[i_axis][i_point].first) + ")");

      if (m_points[i_axis][i_point+1].first <= m_points[i_axis][i_point].first)
        throw std::runtime_error("EtaDistribution: biased positions for axis " + std::to_string(i_axis) + " are not monotonically increasing for points " + std::to_string(i_point) + " (" + std::to_string(m_points[i_axis][i_point].first) + ", " + std::to_string(m_points[i_axis][i_point].second) + ") and " + std::to_string(i_point+1) + " (" + std::to_string(m_points[i_axis][i_point+1].first) + ", " + std::to_string(m_points[i_axis][i_point+1].second) + ")");
    }


    // if there is no point at biases pos = 0, extrapolate one
    // this assumes the function wraps around the pixel edge (which is true for the eta function)
    if (m_points[i_axis].front().first > 1.e-6f) {
      const float x0 = m_points[i_axis].back().first-1;
      const float y0 = m_points[i_axis].back().second-1;
      const float x1 = m_points[i_axis].front().first;
      const float y1 = m_points[i_axis].front().second;
      const float slope = (y1 - y0) / (x1 - x0);
      const float y = y0 - slope * x0; // linear extrapolation to x=0
      m_points[i_axis].insert(m_points[i_axis].begin(), {0.f, y});
    }

    // pre-compute the slope (dTruthPos/dCoG) between each pair of points
    for (size_t i_point=0; i_point < m_points[i_axis].size()-1; ++i_point) {
      const auto& [x0, y0] = m_points[i_axis][i_point];
      const auto& [x1, y1] = m_points[i_axis][i_point+1];
      // m_slopes[i_axis].push_back(y1 > y0 ? (x1 - x0) / (y1 - y0) : 0.f); // 0 for degenerate (flat) segments
      m_slopes[i_axis].push_back( (y1 - y0) / (x1 - x0) );
    }
    // last point needs a slope too. Get this by assuming the function wraps, like above
    auto [x0, y0] = m_points[i_axis].back();
    const auto& [x1, y1] = m_points[i_axis].front();
    x0 += 1.f; y0 += 1.f; // use wrapping-assumption
    // m_slopes[i_axis].push_back(y1 > y0 ? (x1 - x0) / (y1 - y0) : 0.f);
    m_slopes[i_axis].push_back( (y1 - y0) / (x1 - x0) );
  } // end loop axes (u/v)
}

float EtaDistribution::CorrectPos(const int axis, const float biasedPos) const {
  if (axis != 0 && axis != 1)
    throw std::runtime_error("VTXdigi_tools::EtaDistribution::Interpolate: axis must be 0 or 1");
  if (biasedPos < 0.f || biasedPos > 1.f)
    throw std::runtime_error("VTXdigi_tools::EtaDistribution::Interpolate: biasedPos must be in [0, 1] relative to centre of left pixel");

  // find the bracket for biasedPos among the point's x-values
  size_t i_point = 0; // biasedPos lies in bracket (i_point, i_point+1)
  for (size_t i=1; i < m_points.at(axis).size(); ++i) {
    if (m_points.at(axis)[i].first > biasedPos)
      break;
    i_point = i;
  }

  const float slope = m_slopes.at(axis).at(i_point);
  const float correctedPos = m_points.at(axis).at(i_point).second + slope * (biasedPos - m_points.at(axis).at(i_point).first);
  return correctedPos;
}

const std::pair<float, float>& EtaDistribution::GetFunctionBinValue(const int axis, const unsigned int bin) const {
  if (axis != 0 && axis != 1)
    throw std::runtime_error("VTXdigi_tools::EtaDistribution::GetFunctionBinValue: axis must be 0 or 1");
  return m_points.at(axis).at(bin);
}

/* -- Clusterization -- */

std::array<float, 2> Cluster::ComputeCoG() const {
  std::array<float, 2> pos{0.f, 0.f};
  for (const Pixel* pix : pixels) {
    pos[0] += pix->index[0] * pix->charge;
    pos[1] += pix->index[1] * pix->charge;
  }
  pos[0] /= charge;
  pos[1] /= charge;
  return pos;
}

std::array<float, 2> Cluster::ComputeCoG_EtaCorrected(const EtaDistribution& etaDistrib) const {
  std::array<float, 2> pos{0.f, 0.f};

  for (int axis=0; axis<2; ++axis) {
    const int clstLength = GetSize(axis);
    if (clstLength == 1) {
      // no eta correction
      pos[axis] = pixels.front()->index[axis];
    }
    else {
      // apply eta correction for clusters of length >= 2 (for longer clusters, only consider the first and last pixels along this axis, a la https://cds.cern.ch/record/687475/files/note02_049.pdf)

      // find pixels in first and last bin along this axis
      int i_first=std::numeric_limits<int>::max(), i_last=std::numeric_limits<int>::min();
      for (const Pixel* pix : pixels) {
        i_first = std::min(i_first, pix->index[axis]);
        i_last = std::max(i_last, pix->index[axis]);
      }

      // sum up charge in first and last bin
      float charge_first=0.f, charge_last=0.f;
      for (const Pixel* pix : pixels) {
        if (pix->index[axis] == i_first)
          charge_first += pix->charge;
        else if (pix->index[axis] == i_last)
          charge_last += pix->charge;
      }

      // compute pos offset along the axis (as if the cluster was 2 pixels long)
      const float biasedOffset = charge_last / (charge_first + charge_last); // charge center of gracity offset in [0, 1], relative to the center of the first pixel
      const float correctedOffset = etaDistrib.CorrectPos(axis, biasedOffset); // apply eta correction

      pos[axis] = (static_cast<float>(i_first + i_last)) / 2.f - 0.5f + correctedOffset; //
    }
  }
  return pos;
}

int Cluster::GetSize(const int axis) const {
  int min = std::numeric_limits<int>::max();
  int max = std::numeric_limits<int>::min();

  if (axis == 0) { // u
    for (const Pixel* pix : pixels) {
      min = std::min(min, pix->index[0]);
      max = std::max(max, pix->index[0]);
    }
  }
  else if (axis == 1) { // v
    for (const Pixel* pix : pixels) {
      min = std::min(min, pix->index[1]);
      max = std::max(max, pix->index[1]);
    }
  }
  else {
    throw std::runtime_error("Cluster::GetClusterSize: axis must be 0 (u) or 1 (v), got " + std::to_string(axis));
  }

  return max - min + 1; // +1 because of counting: if min=max, cluster size is 1, not 0
}

std::array<float, 2> Cluster::ComputeCoGUncertainty(const std::array<float, 2>& clusterPos) const {
  float sig2_u=0.f, sig2_v=0.f;
  for (const Pixel* pix : pixels) {
    float du = (pix->index[0] - clusterPos[0]);
    float dv = (pix->index[1] - clusterPos[1]);
    sig2_u += pix->charge * du * du;
    sig2_v += pix->charge * dv * dv;
  }
  sig2_u /= charge;
  sig2_v /= charge;
  return {std::sqrt(sig2_u), std::sqrt(sig2_v)};
}

float Cluster::GetSeedPixelCharge() const {
  float maxCharge = 0;
  for (const auto pixel : pixels) {
    if (pixel->charge > maxCharge) {
      maxCharge = pixel->charge;
    }
  }
  return maxCharge;
}


std::array<std::array<int, 2>, 4> GetDirectNeighbors(const std::array<int, 2>& i_uv) {
  return {{
    {i_uv[0] - 1, i_uv[1]}, // left
    {i_uv[0] + 1, i_uv[1]}, // right
    {i_uv[0], i_uv[1] - 1}, // down
    {i_uv[0], i_uv[1] + 1}  // up
  }};
}
std::array<std::array<int, 2>, 8> GetNeighbors(const std::array<int, 2>& i_uv) {
  return {{
    {i_uv[0] - 1, i_uv[1]}, // left
    {i_uv[0] - 1, i_uv[1] + 1}, // upper left
    {i_uv[0], i_uv[1] + 1}, // up
    {i_uv[0] + 1, i_uv[1] + 1}, // upper right
    {i_uv[0] + 1, i_uv[1]}, // right
    {i_uv[0] + 1, i_uv[1] - 1}, // lower right
    {i_uv[0], i_uv[1] - 1}, // lower
    {i_uv[0] - 1, i_uv[1] - 1}, // lower left
  }};
}


std::vector<Cluster> HitMap::ComputeClusters_singePixels() const {
  std::vector<Cluster> clusters;

  for (const auto& p : m_pixels) {
    const Pixel* pixel = &(p.second); // cluster stores pointers to pixels (pixels are stored in HitMap as values)

    clusters.emplace_back();
    clusters.back().pixels.push_back(pixel);
    clusters.back().charge = pixel->charge;
    for (const SimHitWrapper* simHitWrapper : pixel->simHits) {
      clusters.back().simHits.insert(simHitWrapper);
    }
  } // loop over pixelHits

  return clusters;
}

std::vector<Cluster> HitMap::ComputeClusters() const {
  /* Breadth First Search (BFS) implementation for clustering */

  std::vector<Cluster> clusters;
  std::unordered_set<std::array<int, 2>, Hash_PixIndex> visited;

  for (const auto& p : m_pixels) {
    const std::array<int, 2> seed_uv = p.first;
    if (visited.contains(seed_uv))
      continue;

    clusters.emplace_back(); // create new cluster
    clusters.back().pixels.reserve(10); // 10 should include >90% of clusters. i guess.

    std::queue<std::array<int, 2>> queue;
    queue.push(seed_uv);
    visited.insert(seed_uv);

    while (!queue.empty()) {
      const std::array<int, 2> current_uv = queue.front();
      queue.pop();

      /* Add pixl to cluster */
      const Pixel* pixel = &(m_pixels.at(current_uv)); // get pixel pointer from map
      clusters.back().pixels.push_back(pixel);
      clusters.back().charge += pixel->charge;
      for (const SimHitWrapper* simHitWrapper : pixel->simHits) {
        clusters.back().simHits.insert(simHitWrapper);
      }
      /* Add all neighboring pixels to queue */
      // for (const auto& neighbor_uv : GetDirectNeighbors(current_uv)) {
      for (const auto& neighbor_uv : GetNeighbors(current_uv)) {
        if (!m_pixels.contains(neighbor_uv))
          continue;
        if (visited.contains(neighbor_uv))
          continue;
        queue.push(neighbor_uv);
        visited.insert(neighbor_uv);
      }

    } // loop over queue
  } // loop over cluster-seeds
  return clusters;
}

} // namespace VTXdigi_tools
