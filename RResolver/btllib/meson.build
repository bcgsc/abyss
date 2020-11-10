project('btllib', 'cpp',
        version : '1.0.5',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++11' ])

global_args = [ '-pthread' ]
global_link_args = [ '-ldl', '-lpthread' ]
if build_machine.system() == 'darwin'
  global_args += [ '-Xpreprocessor', '-fopenmp' ]
  global_link_args += [ '-Xpreprocessor', '-fopenmp', '-lomp' ]
else
  global_args += [ '-fopenmp' ]
  global_link_args += [ '-fopenmp' ]
endif
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

clang_format = find_program('clang-format', required: false)
if clang_format.found()
  run_target('format', command: 'scripts/format')
endif

swig = find_program('swig', required: false)
if swig.found()
  run_target('wrap', command: 'scripts/wrap')
endif

clang_tidy = find_program('clang-tidy', required: false)
if clang_tidy.found()
  run_target('tidycheck', command: 'scripts/tidycheck')
endif

diff = find_program('diff', required : false)
if diff.found() and clang_format.found()
  run_target('formatcheck', command: 'scripts/formatcheck')
endif

cppcheck = find_program('cppcheck', required: false)
if cppcheck.found()
  run_target('cppcheck', command: 'scripts/cppcheck')
endif

doxygen = find_program('doxygen', required: false)
if doxygen.found()
  run_target('docs', command: 'scripts/docs')
endif

run_target('complete', command: 'scripts/complete')

subdir('include')
subdir('tests')
subdir('extras')
