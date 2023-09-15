src_files = [
    'src/zone.C',
    'src/quantity.C',
    'src/element.C',
    'src/flowequationset.C',
    'src/geometryreference.C',
    'src/family.C',
    'src/familybc.C',
    'src/flowsolution.C',
    'src/base.C',
    'src/range.C',
    'src/dataarray.C',
    'src/bcproperty.C',
    'src/bc.C',
    'src/zonegridconnectivity.C',
    'src/userdefineddata.C',
    'src/gridconnectivity.C',
    'src/bcdataset.C',
    'src/gridcoordinates.C',
    'src/gridconnectivityproperty.C',
    'src/structure_t.C',
    'src/checkerror.C',
    'src/node.C',
    'src/zonebc.C',
    'src/file.C'
]

# Find all header files ending with '.h'
cgnsoo3_H = run_command('find', 'src', '-name', '*.H').stdout().strip().split('\n')

# Installation instructions
lib = library('cgnsoo3', 
               src_files,
               install: true,
               dependencies: cgns_dep)
install_headers(cgnsoo3_H, subdir: 'libcgnsoo3')