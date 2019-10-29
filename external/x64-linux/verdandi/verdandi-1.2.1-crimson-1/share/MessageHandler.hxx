// Copyright (C) 2009, INRIA
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


#ifndef VERDANDI_FILE_SHARE_MESSAGEHANDLER_HXX


namespace Verdandi
{


    //! This class enables objects to communicate using messages.
    class MessageHandler
    {
    public:
        //! Function pointer.
        typedef void (* function_pointer)(void*, string);
        //! Recipient list.
        typedef list<pair<void*, function_pointer> > recipient_list;
        //! Recipient map.
        typedef map<string, recipient_list> recipient_map;

    protected:
        //! Static recipient map.
        static recipient_map recipient_map_;

    public:
        template <class R>
        static void AddRecipient(string recipient, R& object,
                                 function_pointer pointer);
        static void AddRecipient(string recipient, void* object,
                                 function_pointer pointer);
        static void Send(string recipient, string message);
        template <class Sender>
        static void Send(const Sender& sender, string recipient,
                         string message);

        static string GetName();
        template <class R>
        static void RemoveRecipient(R& object);

    private:
        static void SendToList(recipient_list& my_list, string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_MESSAGEHANDLER_HXX
#endif
