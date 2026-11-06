#pragma once

#include <stdexcept>
#include <string>

namespace thermorna::rnastructure {

/**
 * @brief Thrown when parameter files are missing, malformed, or cannot be combined.
 *
 * The message names the file (and line, where there is one) so it can be shown to users as is.
 */
class ParameterError : public std::runtime_error {
   public:
    explicit ParameterError(const std::string& message) : std::runtime_error(message) {}
};

}  // namespace thermorna::rnastructure
