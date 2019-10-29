// Copyright (C) 2008, INRIA
// Author(s): Anne Tilloy
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


#ifndef VERDANDI_FILE_SHARE_LOCKFILE_CXX


#include "LockFile.hxx"

#include <fcntl.h>
#include <cerrno>
#include <cstdio>

namespace Verdandi
{


    //! Creates a lock file.
    /*! If the lock file already exists, this function waits first for this
      file to be removed. If the file is still there after one million
      seconds, the lock fails and this function returns 'false'.
      \param[in] filename name of the lock file.
      \return True if the lock file was successfully created, false otherwise.
    */
    bool Lock(string const& filename)
    {
        int wait = 0;
        const size_t max_wait = 1000000;
#ifdef WIN32
        int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_EXCL);
#else
        int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0200);
#endif
        while (fd < 0)
            if (errno == EEXIST && wait++ < int(max_wait))
            {
#ifdef WIN32
                Sleep(1000000);
                fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_EXCL);
#else
                usleep(1000000);
                fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_EXCL,
                          0200);
#endif
            }
            else
                return false;
        close(fd);
        return true;
    }


    //! Removes a lock file.
    /*!
      \param[in] filename name of the lock file.
    */
    bool Unlock(string const& filename)
    {
        return unlink(filename.c_str()) == 0;
    }


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_LOCKFILE_CXX
#endif
