// VTXdigi_Modular/include/IChargeCollector.h
#pragma once

#include <array>
#include <vector>
#include <optional>

#include "TGeoMatrix.h"

struct VTXdigi_Modular;

namespace VTXdigi_tools {

  class SimHitWrapper; // forward-declare things in include/VTXdigi_tools.h
  class HitMap;

class IChargeCollector {
public:
  virtual ~IChargeCollector() = default;
  virtual void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const = 0;

  /** @brief Compute the eta distribution values from the charge collector
   * @return An optional array containing the eta distribution values in u and v, or std::nullopt if not applicable in the selected implementation
   * @note For an even binning across a single pixel, the vector contains the charge collection centre of gravity for each bin
   * @note Implemeted such to stop LUT implementation details from leaking into the digitizer */
  virtual std::optional< std::array<std::vector<float>, 2> > ComputeEtaDistribution() const { return std::nullopt; }

  float GetChargeCollectionDepthCenter() const { return m_chargeCollectionDepthCenter; }

protected:
  explicit IChargeCollector(const VTXdigi_Modular& digitizer) : m_digitizer(digitizer) {}

  const VTXdigi_Modular& m_digitizer;
  float m_chargeCollectionDepthCenter=0; // Defines the vertical center of the charge collection region in the sensitive volume. Needed for correct digiHit positions (and residual plots) with LUTs where the charge collection varies along their depth (like TPSCo 65nm CIS). in mm, wrt. the sensor local w coordinate (w=0 at center of sensitive volume)
};

std::unique_ptr<IChargeCollector> CreateChargeCollector(const VTXdigi_Modular& digitizer, const std::string& algorithm);

} // namespace VTXdigi_tools
