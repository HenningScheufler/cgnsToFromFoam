project('libcgnsoo_3.0', 'c', 'cpp',
        default_options: ['buildtype=release',
                          'warning_level=0',
                          'c_std=c99',
                          'cpp_std=c++11'])

# mpi flag
mpi = false

if get_option('custom-mpi')
  mpi_dep = []
  mpi = true
# Otherwise they are automatically determined
else
  mpi_dep = [dependency('mpi', language:'c', required : get_option('with-mpi')),
             dependency('mpi', language:'cpp', required : get_option('with-mpi'))]
  if mpi_dep[0].found() or mpi_dep[1].found()
    mpi = true
  endif
endif


subdir('cgns')
subdir('libcgnsoo_3.0')

