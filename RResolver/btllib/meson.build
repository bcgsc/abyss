project('btllib', 'cpp',
        version : '1.0.6',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++11', 'warning_level=3', 'werror=true' ])

global_args = [ '-pthread' ]
global_link_args = [ '-ldl', '-pthread', '-lpthread' ]
if build_machine.system() == 'darwin'
  global_args += [ '-Xpreprocessor', '-fopenmp' ]
  global_link_args += [ '-Xpreprocessor', '-fopenmp', '-lomp' ]
else
  global_args += [ '-fopenmp' ]
  global_link_args += [ '-fopenmp' ]
endif

sdsl = meson.get_compiler('cpp').find_library('sdsl', required : false, dirs : meson.source_root()+ '/external/sdsl-lite/installdir/lib')
global_args += [ '-I' + meson.source_root() + '/external/sdsl-lite/installdir/include' ]
global_link_args += [ '-L' + meson.source_root() + '/external/sdsl-lite/installdir/lib', '-lsdsl', '-ldivsufsort', '-ldivsufsort64' ]

add_global_arguments(global_args, language : 'cpp')
add_global_link_arguments(global_link_args, language : 'cpp')

compiler = meson.get_compiler('cpp')
compiler_id = compiler.get_id()
compiler_version = compiler.version()
if compiler_id == 'gcc' and compiler_version.version_compare('<4.8.1')
  error('GCC ' + compiler_version + ' doesn\'t fully support C++11')
endif
if compiler_id == 'clang' and compiler_version.version_compare('<3.3.0')
  error('Clang ' + compiler_version + ' doesn\'t fully support C++11')
endif

btllib_include = include_directories('include')

rootpath = meson.source_root()

run_target('build-sdsl', command: join_paths(rootpath, 'scripts/build_sdsl'))

clang_format = find_program('clang-format', required : false)
if clang_format.found()
  run_target('format', command: join_paths(rootpath, 'scripts/format'))
endif

swig = find_program('swig', required : false)
if swig.found()
  run_target('wrap', command: join_paths(rootpath, 'scripts/wrap'))
endif

clang_tidy = find_program('clang-tidy', required : false)
if clang_tidy.found()
  run_target('tidycheck', command: [ join_paths(rootpath, 'scripts/tidycheck') ] + global_args + global_link_args)
endif

diff = find_program('diff', required : false)
if diff.found() and clang_format.found()
  run_target('formatcheck', command: join_paths(rootpath, 'scripts/formatcheck'))
endif

cppcheck = find_program('cppcheck', required : false)
if cppcheck.found()
  run_target('cppcheck', command: join_paths(rootpath, 'scripts/cppcheck'))
endif

doxygen = find_program('doxygen', required : false)
if doxygen.found()
  run_target('docs', command: join_paths(rootpath, 'scripts/docs'))
endif

run_target('complete', command: join_paths(rootpath, 'scripts/complete'))

subdir('include')
subdir('tests')
subdir('extras')
