pyflags_cmd = run_command('../../scripts/get_python_flags')

if pyflags_cmd.stderr() != ''
  error(pyflags_cmd.stderr())
endif

pyflags = pyflags_cmd.stdout().split('\n')

if pyflags[0] == 'ERROR'
  error(pyflags[1].strip())
endif

cpp_extra_args = pyflags[1].strip().split()
link_extra_args = pyflags[2].strip().split()

cpp_extra_args += ['-w']

# Ensure proper ordering, arguments will be deduplicated anyway
cpp_extra_args = global_args + cpp_extra_args
link_extra_args = global_link_args + link_extra_args

py_sources = ['btllib_wrap.cxx']
btllib_py = shared_module('_btllib',
                           py_sources,
                           name_suffix : 'so',
                           include_directories : btllib_include,
                           dependencies : deps,
                           cpp_args : cpp_extra_args,
                           link_args : link_extra_args,
                           install : true,
                           install_dir : 'python',
                           name_prefix : '')

python_files = ['btllib.py']
install_data(python_files, install_dir : 'python')