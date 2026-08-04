/// SPDX-License-Identifier: GPL-3.0-or-later
/// SPDX-FileCopyrightText: 2026 Frederik Hennig <frederik.hennig@fau.de>

#pragma once

#include <stdexcept>

namespace walberla::sweepgen
{

/**
 * @brief Designates a violated constraint to a sweep parameter
 */
class SweepConstraintError : public std::invalid_argument
{
   using std::invalid_argument::invalid_argument;
};

} // namespace walberla::sweepgen