#ifndef PAPPUS_CONTRACTS_HPP
#define PAPPUS_CONTRACTS_HPP

#include <libassert/assert.hpp>

// Matches Operon's operon/core/contracts.hpp macro text exactly, so a TU
// that includes both never sees a macro-redefinition warning.
#define ENSURE ASSERT
#define EXPECT ASSERT

#endif
