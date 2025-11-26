#include "ChargeCollector.h"

namespace VTXdigi_details {

void ChargeCollector_LUT::Collect() const {

}

void ChargeCollector_Drift::Collect() const {

}

std::unique_ptr<IChargeCollector> CreateChargeCollector(const std::string& algorithm) {
  if (algorithm == "LookupTable_Allpix2") {
    return std::make_unique<ChargeCollector_LUT>();
  } else if (algorithm == "Drift") {
    return std::make_unique<ChargeCollector_Drift>();
  } else {
    throw std::runtime_error("Unknown ChargeCollector type: " + algorithm);
  }
}

} // namespace VTXdigi_details


