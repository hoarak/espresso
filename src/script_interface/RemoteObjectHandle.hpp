/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SCRIPT_INTERFACE_PARALLEL_SCRIPT_INTERFACE_SLAVE_HPP
#define SCRIPT_INTERFACE_PARALLEL_SCRIPT_INTERFACE_SLAVE_HPP

#include "MpiCallbacks.hpp"
#include "ObjectHandle.hpp"

#include "ParallelScriptInterface.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/utility.hpp>
#include <jmorecfg.h>

namespace ScriptInterface {
class RemoteObjectHandle {
private:
  const auto &comm() const { return m_callback_id.cb()->comm(); }

  using CallbackAction = ParallelScriptInterface::CallbackAction;

public:
  RemoteObjectHandle(Communication::MpiCallbacks *cb)
      : m_callback_id(cb, [this](CallbackAction a) { mpi_slave(a); }) {}

  std::shared_ptr<ObjectHandle> m_p;

  static auto &get_translation_table() {
    static std::map<ObjectId, std::weak_ptr<ObjectHandle>> m_translation_table;

    return m_translation_table;
  }

  /* If the variant encapsulates an object id we translate the
   master id to a local one */
  static void translate_id(Variant &v) {
    if (is_type<ObjectId>(v)) {
      v = get_translation_table().at(get_value<ObjectId>(v)).lock()->id();
    }
  }

  static void translate_id(VariantMap &map) {
    for (auto &e : map) {
      translate_id(e.second);
    }
  }

  VariantMap bcast_variant_map() const {
    VariantMap ret;
    boost::mpi::broadcast(comm(), ret, 0);

    /* If the parameter is a object we have to translate it first to a
       local id.
    */
    translate_id(ret);

    return ret;
  }

private:
  Communication::CallbackHandle<CallbackAction> m_callback_id;
  void mpi_slave(CallbackAction action) {
    switch (action) {
    case CallbackAction::CONSTRUCT: {
      std::pair<ObjectId, std::string> what;
      boost::mpi::broadcast(comm(), what, 0);
      auto const parameters = bcast_variant_map();

      m_p = ObjectHandle::make_shared(
          what.second, ObjectHandle::CreationPolicy::LOCAL, parameters);

      get_translation_table()[what.first] = m_p;
      break;
    }
    case CallbackAction::SET_PARAMETER: {
      std::pair<std::string, Variant> d;
      boost::mpi::broadcast(comm(), d, 0);

      /* If the parameter is a object we have to translate it first to a
         local id.
      */
      translate_id(d.second);
      m_p->set_parameter(d.first, d.second);
      break;
    }
    case CallbackAction::CALL_METHOD: {
      /* Name of the method and parameters */
      std::pair<std::string, VariantMap> d;

      /* Broadcast method name and parameters */
      boost::mpi::broadcast(comm(), d, 0);

      translate_id(d.second);

      /* Forward to the local instance. */
      m_p->call_method(d.first, d.second);

      break;
    }
    case CallbackAction::DELETE: {
      delete this;
      break;
    }
    }
  }
};
} /* namespace ScriptInterface */

#endif
