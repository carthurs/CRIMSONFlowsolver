import distutils.sysconfig
env = Environment(SWIGFLAGS = ['-Wall', '-c++', '-python'],
                  CPPPATH = [distutils.sysconfig.get_python_inc()],
                  SHLIBPREFIX = "")
env.Append(CPPFLAGS = "-DSELDON_DEBUG_LEVEL_4")
env.SharedLibrary('_seldon.so', ['Seldon.cpp', 'seldon.i'])
