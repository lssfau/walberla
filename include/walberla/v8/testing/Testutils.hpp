//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Testutils.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <ranges>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <string_view>

#include "walberla/v8/Memory.hpp"
#include "walberla/v8/Sweep.hpp"

namespace walberla::v8::testing
{

/**
 * @defgroup v8core-testutils Testsuite Utilities (v8::testing)
 * @brief Unit testing toolkit for waLBerla
 * 
 * This toolkit of test utilities is designed to facilitate unit testing
 * of the waLBerla framework.
 * 
 * See also [Testing in the V8 core contributors guide](#v8-contrib-testing).
 * 
 * ## Creating a Test Executable
 * 
 * Tests are grouped into test executables. We thematically group as many tests
 * as sensible into a single executable to accelerate compilation times.
 * Test executables follow the naming scheme `TestX.cpp`.
 * To register a test executable with CMake, use the `waLBerla_add_test_executable`
 * function:
 * 
 * ```
 * waLBerla_add_test_executable( TestX TestX.cpp )
 * ```
 * 
 * ## Registering Tests
 * 
 * Every test executable must include the following scaffolding, where the tests are registered
 * with a `TestsRunner` instance:
 * 
 * ```
 * #include "walberla/v8/Testing.hpp"
 * 
 * int main(int argc, char** argv) {
 *    walberla::mpi::Environment env{ argc, argv };
 * 
 *    return walberla::v8::testing::TestsRunner({ 
 *         // TEST FUNCTIONS
 *         {"MyTest1", &myTest1},
 *         {"MyTest2", &myTest2},
 *         // ...
 *    }).run(argc, argv);
 * }
 * ```
 * 
 * The individual tests (`MyTest1`, `MyTest2`, ...) are zero-argument `void` functions.
 * To automatically run the tests, they must further be registered with `CTest` in the `CMakeLists.txt`,
 * using `walberla_v8_add_tests`:
 * 
 * ```
 * walberla_v8_add_tests( TestX IDS MyTest1 MyTest2 ... )
 * ```
 * 
 * Here, `TestX` is the name of the test application, the `IDS` are the test names as passed to the `TestsRunner`.
 * 
 * ## Writing Tests
 * 
 * Tests should be written using the assertion functions from this module, listed below.
 * 
 */

namespace printing
{
std::string format_range(const std::ranges::range auto& range)
{
   std::stringstream str;
   str << "[";
   bool first = true;
   for (auto x : range)
   {
      if (!first) { str << ", "; }
      else
      {
         first = false;
      }
      str << std::format("{}", x);
   }
   str << "]";
   return str.str();
}

std::string format_range(const auto& not_a_range)
{
   std::stringstream str;
   str << not_a_range;
   return str.str();
}

} // namespace printing

class AssertionError : public std::exception
{
 public:
   AssertionError(const std::string& msg, const std::source_location& loc) : msg_{ msg }, loc_{ loc } {}

   [[nodiscard]] const char* what() const noexcept override { return msg_.c_str(); }

   const std::string& message() const { return msg_; }

   const std::source_location& location() const { return loc_; }

 private:
   std::string msg_;
   std::source_location loc_;
};

/**
 * @brief Primary test orchestrator
 * @ingroup v8core-testutils
 */
class TestsRunner
{
 public:
   using TestFunction = std::function< void() >;

   /**
    * @brief Create a test runner instance, registering a set of test functions
    */
   TestsRunner(std::initializer_list< std::tuple< std::string, TestFunction > > init_list)
   {
      for (const auto& [k, v] : init_list)
      {
         tests_[k] = v;
      }
   }

   /**
    * @brief Run a single test by its ID
    */
   int runTest(const std::string& id)
   {
      if (auto func = tests_.find(id); func != tests_.end())
      {
         try
         {
            (std::get< 1 >(*func))();
            return EXIT_SUCCESS;
         } catch (const AssertionError& err)
         {
            auto& loc = err.location();

            std::stringstream errStream;
            errStream << "Assertion failed: " << err.message() << "\n"
                      << "    At " << loc.file_name() << '(' << loc.line() << ':' << loc.column() << ")\n"
                      << "    in function `" << loc.function_name() << "`" << "\n";
            auto logMessage = errStream.str();

            walberla::Abort::instance()->abort(logMessage, err.location().file_name(), (int) err.location().line());
            return 1;
         }
      }
      else
      {
         std::clog << "Invalid test ID: " << id << "\n";
         return 2;
      }
   }

   /**
    * @brief Run tests as selected from the command line arguments
    */
   int run(int argc, char** argv)
   {
      const std::span< char* > args(argv, size_t(argc));

      if (args.size() < 2)
      {
         std::clog << "No test ID was provided\n";
         std::exit(2);
      }

      const std::string testId{ args[1] };
      return runTest(testId);
   }

 private:
   std::map< std::string, TestFunction > tests_;
};

/**
 * @brief Create and obtain a unique temporary directory.
 * @ingroup v8core-testutils
 * 
 * @note To override the location where temporary directories are placed,
 *       set the `TMPDIR` environment variable.
 *       See also [`std::filesystem::temp_directory_path`](https://en.cppreference.com/w/cpp/filesystem/temp_directory_path.html)
 */
inline std::filesystem::path tmp_dir()
{
   const auto now     = std::chrono::steady_clock::now();
   const auto dirname = std::format("walberla-testing-{}", now.time_since_epoch().count());
   const auto tmpDir  = std::filesystem::temp_directory_path() / dirname;
   std::filesystem::create_directories(tmpDir);
   return tmpDir;
}

/**
 * @brief Check if the given boolean condition is `true`.
 * @ingroup v8core-testutils
 */
inline void assert_true(bool cond, const std::source_location loc = std::source_location::current())
{
   if (!cond) { throw AssertionError{ std::format("Assertion failed: Condition was false, expected true"), loc }; }
}

/**
 * @brief Check if the given boolean condition is `false`.
 * @ingroup v8core-testutils
 */
inline void assert_false(bool cond, const std::source_location loc = std::source_location::current())
{
   if (cond) { throw AssertionError{ std::format("Assertion failed: Condition was true, expected false"), loc }; }
}

/**
 * @brief Check if the given callback function raises an expected exception type.
 * @ingroup v8core-testutils
 * 
 * Example:
 * 
 * ```
 * throws< std::invalid_argument >([&](){ 
 *    // code that should throw ... 
 * });
 * ```
 */
template< typename TException >
void throws(const std::function< void() >& func, const std::source_location loc = std::source_location::current())
{
   bool caught{ false };

   //  NOLINTBEGIN(bugprone-empty-catch)
   try
   {
      func();
   } catch (const TException&) { caught = true; } catch (const std::exception&)
   {}
   //  NOLINTEND(bugprone-empty-catch)

   if (!caught) { throw AssertionError("Expected an exception, but none (or the wrong type) was thrown", loc); }
}

/**
 * @brief Check if two values are equal.
 * @ingroup v8core-testutils
 * @note For floating-point comparisons, use `assert_close` instead.
 */
template< std::equality_comparable T >
void assert_equal(const T& actual, const T& desired, const std::source_location loc = std::source_location::current())
{
   if (actual != desired)
   {
      throw AssertionError(std::format("Values were not equal. Actual: {}, Desired: {}", printing::format_range(actual),
                                       printing::format_range(desired)),
                           loc);
   }
}

/**
 * @brief Check if two values are inequal.
 * @ingroup v8core-testutils
 */
template< std::equality_comparable T >
void assert_inequal(const T& left, const T& right, const std::source_location loc = std::source_location::current())
{
   if (left == right)
   {
      throw AssertionError(
         std::format("{} was equal to {}", printing::format_range(left), printing::format_range(right)), loc);
   }
}

/**
 * @brief Check if one value is less than another.
 * @ingroup v8core-testutils
 */
template< std::totally_ordered T >
void assert_less(const T& left, const T& right, const std::source_location loc = std::source_location::current())
{
   if (!(left < right))
   {
      throw AssertionError(
         std::format("{} not less than {}", printing::format_range(left), printing::format_range(right)), loc);
   }
}

/**
 * @brief Check if one value is greater than another.
 * @ingroup v8core-testutils
 */
template< std::totally_ordered T >
void assert_greater(const T& left, const T& right, const std::source_location loc = std::source_location::current())
{
   if (!(left > right))
   {
      throw AssertionError(
         std::format("{} not greater than {}", printing::format_range(left), printing::format_range(right)), loc);
   }
}

/**
 * @brief Check if one value is less or equal to another.
 * @ingroup v8core-testutils
 */
template< std::totally_ordered T >
void assert_less_equal(const T& left, const T& right, const std::source_location loc = std::source_location::current())
{
   if (!(left <= right))
   {
      throw AssertionError(
         std::format("{} not less than or equal to {}", printing::format_range(left), printing::format_range(right)),
         loc);
   }
}

/**
 * @brief Check if one value is greater or equal to another.
 * @ingroup v8core-testutils
 */
template< std::totally_ordered T >
void assert_greater_equal(const T& left, const T& right,
                          const std::source_location loc = std::source_location::current())
{
   if (!(left >= right))
   {
      throw AssertionError(
         std::format("{} not greater than or equal to {}", printing::format_range(left), printing::format_range(right)),
         loc);
   }
}

/**
 * @brief Check if two ranges have identical elements
 * @ingroup v8core-testutils
 */
template< std::ranges::range R1, std::ranges::range R2 >
   requires(std::equality_comparable_with< std::ranges::range_value_t< R1 >, std::ranges::range_value_t< R2 > >)
void assert_range_equal(const R1& actual, const R2& desired,
                        const std::source_location loc = std::source_location::current())
{
   if (!std::ranges::equal(actual, desired))
   {
      throw AssertionError(std::format("Ranges were not equal. Actual: {}, Desired: {}", printing::format_range(actual),
                                       printing::format_range(desired)),
                           loc);
   }
}

/**
 * @brief Check if all elements of a range are equal to a desired value
 * @ingroup v8core-testutils
 * @note If the range is empty, this assertion will succeed.
 */
template< std::ranges::range R1 >
void assert_range_equal(const R1& actual, const std::ranges::range_value_t< R1 >& desired,
                        const std::source_location loc = std::source_location::current())
{
   if (std::ranges::any_of(actual, [desired](auto& v) { return v != desired; }))
   {
      throw AssertionError(std::format("Ranges were not equal. Actual: {}, Desired: {}", printing::format_range(actual),
                                       printing::format_range(desired)),
                           loc);
   }
}

template< std::ranges::range R1, std::ranges::range R2 >
   requires(std::equality_comparable_with< std::ranges::range_value_t< R1 >, std::ranges::range_value_t< R2 > >)
void assert_range_inequal(R1 actual, R2 desired, const std::source_location loc = std::source_location::current())
{
   if (std::ranges::equal(actual, desired))
   {
      throw AssertionError(
         std::format("{} was equal to {}", printing::format_range(actual), printing::format_range(desired)), loc);
   }
}

/**
 * @brief Floating-point near-equality assertions with given tolerances
 * @ingroup v8core-testutils
 * 
 * For given `atol` and `rtol` values, `is_close(actual, desired)` succeeds according to:
 * ```
 * tolerance = atol + rtol * abs(desired);
 * success = abs(actual - desired) <= tolerance;
 * ```
 * 
 * Example:
 * ```
 * with_tolerance(1e-12, 1e-6).assert_close(actual, desired);
 * ```
 */
class with_tolerance
{
 private:
   double atol_;
   double rtol_;

   template< std::floating_point T >
   bool isclose(T actual, T desired) const
   {
      T tolerance{ T(atol_) + T(rtol_) * std::abs(desired) };
      return std::abs(actual - desired) <= tolerance;
   }

   [[noreturn]] void fail_not_close(auto& actual, auto& desired, const std::source_location loc) const
   {
      throw AssertionError(std::format("Values not equal to atol={}, rtol={}. Actual: {}, Desired: {}", atol_, rtol_,
                                       printing::format_range(actual), printing::format_range(desired)),
                           loc);
   }

 public:
   /**
    * @brief Set the comparison tolerances
    */
   with_tolerance(double atol, double rtol) : atol_{ atol }, rtol_{ rtol } {}

   /**
    * @brief Check if two values are equal up to the set tolerances
    */
   template< std::floating_point T >
   void assert_close(T actual, T desired, const std::source_location loc = std::source_location::current()) const
   {
      if (!isclose(actual, desired)) { fail_not_close(actual, desired, loc); }
   }

   /**
    * @brief Check if the elements of two ranges are pairwise equal up to the set tolerances.
    * @note Assertion will fail if the ranges are not of the same length.
    */
   template< std::ranges::sized_range R1, std::ranges::sized_range R2 >
      requires(std::floating_point< std::ranges::range_value_t< R1 > > &&
               std::same_as< std::ranges::range_value_t< R1 >, std::ranges::range_value_t< R2 > >)
   void assert_allclose(const R1& actual, const R2& desired,
                        const std::source_location loc = std::source_location::current()) const
   {
      const size_t actualSize{ std::ranges::size(actual) };
      const size_t desiredSize{ std::ranges::size(desired) };

      if (actualSize != desiredSize) { fail_not_close(actual, desired, loc); }

      bool do_fail{ false };

      for (auto [aIt, dIt] = std::make_tuple(std::ranges::begin(actual), std::ranges::begin(desired));
           aIt != std::ranges::end(actual) && dIt != std::ranges::end(desired); ++aIt, ++dIt)
      {
         auto& actualValue  = *aIt;
         auto& desiredValue = *dIt;

         do_fail = do_fail || !this->isclose(actualValue, desiredValue);

         if (do_fail) { break; }
      }

      if (do_fail) { fail_not_close(actual, desired, loc); }
   }

   /**
    * @brief Check if the elements of a range are equal to a desired value up to the set tolerances.
    * @note Assertion will succeed if the range is empty
    */
   template< std::ranges::sized_range R1 >
      requires(std::floating_point< std::ranges::range_value_t< R1 > >)
   void assert_allclose(const R1& actual, std::ranges::range_value_t< R1 > desired,
                        const std::source_location loc = std::source_location::current()) const
   {
      bool do_fail{ false };

      for (auto& actualValue: actual)
      {
         do_fail = do_fail || !this->isclose(actualValue, desired);

         if (do_fail) { break; }
      }

      if (do_fail) { fail_not_close(actual, desired, loc); }
   }

   /**
    * @brief Check if the elements of two waLBerla fields are equal up to the set tolerances.
    */
   template< memory::IFieldView TFieldView >
   void assert_allclose(const TFieldView& actual, const TFieldView& desired, bool withGhostLayers = false,
                        const std::source_location loc = std::source_location::current()) const
   {
      assert_range_equal(actual.shape(), desired.shape(), loc);

      if (withGhostLayers) { assert_equal(actual.numGhostLayers(), desired.numGhostLayers()); }

      const cell_idx_t gls{ withGhostLayers ? cell_idx_c(desired.numGhostLayers()) : 0 };
      const CellInterval ci{ { -gls, -gls, -gls },
                             {
                                cell_idx_c(desired.shape()[0]) + gls - 1, //
                                cell_idx_c(desired.shape()[1]) + gls - 1, //
                                cell_idx_c(desired.shape()[2]) + gls - 1,
                             } };

      const size_t totalEntries{ ci.numCells() * TFieldView::F_SIZE };
      size_t numMismatched{ 0 };

      sweep::forAllCells(exectag::Serial{}, ci, [&](Cell cell) {
         for (cell_idx_t q = 0; q < cell_idx_c(TFieldView::F_SIZE); ++q)
         {
            auto& actualValue{ actual(cell, q) };
            auto& desiredValue{ desired(cell, q) };

            if (!this->isclose(actualValue, desiredValue)) { numMismatched++; }
         }
      });

      if (numMismatched > 0)
      {
         const std::string err =
            std::format("Field entries not equal to atol={}, rtol={}. Number of mismatched entries: {} / {}", atol_,
                        rtol_, numMismatched, totalEntries);
         throw AssertionError(err, loc);
      }
   }
};

/**
 * See `with_tolerance::assert_close`.
 */
template< std::floating_point T >
inline void assert_close(T actual, T desired, const std::source_location loc = std::source_location::current())
{
   with_tolerance(0.0, 1e-7).assert_close(actual, desired, loc);
}

/**
 * See `with_tolerance::assert_allclose`.
 */
template< std::ranges::range R1, std::ranges::range R2 >
   requires(std::floating_point< std::ranges::range_value_t< R1 > > &&
            std::same_as< std::ranges::range_value_t< R1 >, std::ranges::range_value_t< R2 > >)
void assert_allclose(const R1& actual, const R2& desired,
                     const std::source_location loc = std::source_location::current())
{
   with_tolerance(0.0, 1e-7).assert_allclose(actual, desired, loc);
}

/**
 * See `with_tolerance::assert_allclose`.
 */
template< std::ranges::sized_range R1 >
   requires(std::floating_point< std::ranges::range_value_t< R1 > >)
void assert_allclose(const R1& actual, std::ranges::range_value_t< R1 > desired,
                     const std::source_location loc = std::source_location::current())
{
   with_tolerance(0.0, 1e-7).assert_allclose(actual, desired, loc);
}

/**
 * See `with_tolerance::assert_allclose`.
 */
template< memory::IFieldView TFieldView >
void assert_allclose(const TFieldView& actual, const TFieldView& desired, bool withGhostLayers = false,
                     const std::source_location loc = std::source_location::current())
{
   with_tolerance(0.0, 1e-7).assert_allclose(actual, desired, withGhostLayers, loc);
}

} // namespace walberla::v8::testing