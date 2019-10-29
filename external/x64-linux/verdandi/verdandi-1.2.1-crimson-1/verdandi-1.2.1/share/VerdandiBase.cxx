// Copyright (C) 2009 INRIA
// Author(s): Claire Mouton
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_SHARE_VERDANDIBASE_CXX


#include "VerdandiBase.hxx"


namespace Verdandi
{


    //! Default constructor.
    VerdandiBase::VerdandiBase()
    {
    }


    //! Destructor.
    VerdandiBase::~VerdandiBase()
    {
        MessageHandler::RemoveRecipient(*this);
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string VerdandiBase::GetName() const
    {
        return "VerdandiBase";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    void VerdandiBase::Message(string message)
    {
    }


    //! Receives and handles a message with a static method.
    /*
      \param[in] message the received message.
    */
    void VerdandiBase::StaticMessage(void* object, string message)
    {
        reinterpret_cast<VerdandiBase*>(object)->Message(message);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_VERDANDIBASE_CXX
#endif
