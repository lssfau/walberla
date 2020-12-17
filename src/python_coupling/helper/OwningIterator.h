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
//! \file OwningIterator.h
//! \ingroup python_coupling
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include <mutex>
#include <pybind11/pybind11.h>

namespace walberla {
namespace python_coupling {

namespace py = pybind11;

namespace detail {

template <typename T, py::return_value_policy Policy>
struct owning_iterator_state {
   owning_iterator_state(T _obj)
   : obj(_obj), it(obj.begin()), first_or_done(true) {}
   T obj;
   typename T::iterator it;
   bool first_or_done;
   static std::once_flag registered;
};

template <typename T, py::return_value_policy Policy>
std::once_flag owning_iterator_state<T, Policy>::registered;

} // namespace detail

template <py::return_value_policy Policy = py::return_value_policy::reference_internal,
          typename T>
py::iterator make_owning_iterator(T obj) {
   using state = detail::owning_iterator_state<T, Policy>;

   std::call_once(state::registered, []() {
      py::class_<state>(py::handle(), "owning_iterator", py::module_local())
         .def("__iter__", [](state &s) -> state& { return s; })
         .def("__next__", [](state &s) -> typename T::value_type {
            if (!s.first_or_done)
               ++s.it;
            else
               s.first_or_done = false;
            if (s.it == s.obj.end()) {
               s.first_or_done = true;
               throw py::stop_iteration();
            }
            return *s.it;
         }, py::keep_alive< 0, 1 >(), Policy);
   });

   return cast(state(obj));
}

} // namespace python_coupling
} // namespace walberla
