pyflags = run_command('../../scripts/get_python_flags').stdout().split('\n')
cpp_extra_args = pyflags[0].strip().split()
link_extra_args = pyflags[1].strip().split()

cpp_extra_args += ['-w']

py_sources = ['btllib_wrap.cxx']
btllib_py = shared_library('_btllib',
                           py_sources,
                           include_directories : btllib_include,
                           cpp_args : cpp_extra_args,
                           link_args : link_extra_args,
                           install : true,
                           install_dir : 'python',
                           name_prefix : '')

python_files = ['btllib.py']
install_data(python_files, install_dir : 'python')