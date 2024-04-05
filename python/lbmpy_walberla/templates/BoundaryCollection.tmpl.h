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
//! \\file {{class_name}}.h
//! \\author lbmpy
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "domain_decomposition/IBlock.h"

{% if target is equalto 'gpu' -%}
#include "gpu/GPUWrapper.h"
{%- endif %}

{% for include in includes -%}
#include {{include}}
{% endfor %}


namespace walberla{
namespace {{namespace}} {

template <typename FlagField_T>
class {{class_name}}
{
 public:
   enum Type { ALL = 0, INNER = 1, OUTER = 2 };


   {{class_name}}( {{- ["const shared_ptr<StructuredBlockForest> & blocks", "BlockDataID flagID_", "BlockDataID pdfsID_", "FlagUID domainUID_", [kernel_list|generate_constructor_parameters(['indexVector', 'indexVectorSize', 'pdfs'])], additional_constructor_arguments] | type_identifier_list -}} )
      : blocks_(blocks), flagID(flagID_), pdfsID(pdfsID_), domainUID(domainUID_)
   {
      {% for object_name, boundary_class, kernel, additional_data_handler in zip(object_names, boundary_classes, kernel_list, additional_data_handlers) -%}

      {{object_name}} = std::make_shared< {{boundary_class}} >({{- ["blocks", "pdfsID", [kernel|generate_function_collection_call(['indexVector', 'indexVectorSize', 'pdfs', 'timestep', 'gpuStream'], use_field_ids=True)], additional_data_handler.constructor_argument_name] | type_identifier_list -}});
      {% endfor %}

      {% for object_name, flag_uid in zip(object_names, flag_uids) -%}
      {{object_name}}->fillFromFlagField<FlagField_T>(blocks, flagID, walberla::FlagUID("{{flag_uid}}"), domainUID);
      {% endfor %}
   }

   void run ({{- ["IBlock * block", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}})
   {
      {% for object_name in object_names -%}
      {{object_name}}->run({{- ["block", ["stream"] if target == 'gpu' else []] | type_identifier_list -}});
      {% endfor %}
   }

   void inner ({{- ["IBlock * block", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}})
   {
      {% for object_name in object_names -%}
      {{object_name}}->inner({{- ["block", ["stream"] if target == 'gpu' else []] | type_identifier_list -}});
      {% endfor %}
   }

   void outer ({{- ["IBlock * block", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}})
   {
      {% for object_name in object_names -%}
      {{object_name}}->outer({{- ["block", ["stream"] if target == 'gpu' else []] | type_identifier_list -}});
      {% endfor %}
   }

   void operator() ({{- ["IBlock * block", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}})
   {
      run({{- ["block", ["stream"] if target == 'gpu' else []] | type_identifier_list -}});
   }

   std::function<void (IBlock *)> getSweep({{- ["Type type = Type::ALL", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}})
   {
      switch (type)
      {
      case Type::INNER:
         return [{{- ["this", ["stream"] if target == 'gpu' else []] | type_identifier_list -}}](IBlock* block) { this->inner({{- ["block", ["stream"] if target == 'gpu' else []] | type_identifier_list -}}); };
      case Type::OUTER:
         return [{{- ["this", ["stream"] if target == 'gpu' else []] | type_identifier_list -}}](IBlock* block) { this->outer({{- ["block", ["stream"] if target == 'gpu' else []] | type_identifier_list -}}); };
      default:
         return [{{- ["this", ["stream"] if target == 'gpu' else []] | type_identifier_list -}}](IBlock* block) { this->run({{- ["block", ["stream"] if target == 'gpu' else []] | type_identifier_list -}}); };
      }
   }

   weak_ptr< StructuredBlockStorage > blocks_;
   BlockDataID flagID;
   BlockDataID pdfsID;
   walberla::FlagUID domainUID;

   {% for object_name, boundary_class in zip(object_names, boundary_classes) -%}
   shared_ptr<{{boundary_class}}> {{object_name}};
   {% endfor %}
};

}
}
