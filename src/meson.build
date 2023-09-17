project('libcgnsoo_3.0', 'c', 'cpp',
        version: '3.0',
        default_options: ['buildtype=release',
                          'warning_level=0',
                          'c_std=c99',
                          'cpp_std=c++11'])

cgns_dep = dependency('cgns', required : false)

if not cgns_dep.found()
    cgns_proj = subproject('cgns')
    cgns_dep = cgns_proj.get_variable('cgns_dep')
endif

subdir('libcgnsoo_3.0')
