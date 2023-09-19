project('hdf', 'c', 'cpp',
        default_options: ['buildtype=release',
                          'warning_level=0',
                          'c_std=c99',
                          'cpp_std=c++11'])

r = run_command(['cmake',
                '-G Ninja',
                '-S .',
                '-B build',
                '-DCMAKE_INSTALL_PREFIX=' + get_option('prefix'),
                '-DHDF5_ENABLE_PARALLEL=ON',
                '-DHDF5_BUILD_TOOLS=NO',
                '-DHDF5_BUILD_UTILS=No',
                '-DHDF5_BUILD_EXAMPLES=No',
                '-DBUILD_TESTING=No'], check: true)

output = r.stdout().strip()
errortxt = r.stderr().strip()
message('output',output)
message('errortxt',errortxt)

r = run_command(['cmake',
                '--build',
                'build'], check: true)

output = r.stdout().strip()
errortxt = r.stderr().strip()
message('output',output)
message('errortxt',errortxt)

r = run_command(['cmake',
                '--install',
                'build'], check: true)
output = r.stdout().strip()
errortxt = r.stderr().strip()
message('output',output)
message('errortxt',errortxt)