#pragma once

#include <memory>
#include "../include/VTXdigi_Modular.h"

namespace VTXdigi_details {

using ::VTXdigi_Modular;

/* -- Charge collector algorithms -- */

class ChargeCollector_LUT : public IChargeCollector {

  public:
    explicit ChargeCollector_LUT(const VTXdigi_Modular& digitizer);
    void Collect() const override;
};

class ChargeCollector_Drift : public IChargeCollector {

public:
  void Collect() const override;
};

/** @brief Holds position & information about path through the sensor */
struct Path;

Path ComputePath(const dd4hep::rec::ISurface& surface, const edm4hep::SimTrackerHit& simHit);

} // namespace VTXdigi_details