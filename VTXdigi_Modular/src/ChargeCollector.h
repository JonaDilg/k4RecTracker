#pragma once

#include <memory>

namespace VTXdigi_details {

/* -- Charge collector algorithms -- */

class IChargeCollector {
public:
  virtual ~IChargeCollector() = default;

  virtual void Collect() const = 0;
};

class ChargeCollector_LUT : public IChargeCollector {

public:
  void Collect() const override;
};

class ChargeCollector_Drift : public IChargeCollector {

public:
  void Collect() const override;
};

std::unique_ptr<IChargeCollector> CreateChargeCollector(const std::string& algorithm);


} // namespace VTXdigi_details