# Copyright (C) 2010 Vivien Mallet
#
# This file is part of Ops, a library for parsing Lua configuration files.
#
# Ops is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# Ops is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Ops. If not, see http://www.gnu.org/licenses/.


import distutils.sysconfig, os
env = Environment(ENV = os.environ,
                  SWIGFLAGS = ['-Wall', '-c++', '-python'],
                  CPPPATH = [distutils.sysconfig.get_python_inc(),
                             "/usr/include/lua5.1/"],
                  SHLIBPREFIX = "")

conf = Configure(env)
if not conf.CheckLib("lua5.1"):
    conf.CheckLib("lua")

env.Append(CPPFLAGS = "-DOPS_WITH_EXCEPTION")
if env['PLATFORM'] == 'win32':
	env.Append(SHLIBSUFFIX = ".pyd")
	env.Replace(LINK = "LINK")
	env.SharedLibrary('_ops', ['Ops.cpp', 'ops.i'])
else:
    env.SharedLibrary('_ops.so', ['Ops.cpp', 'ops.i'])
