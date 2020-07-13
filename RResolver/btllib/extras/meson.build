python3 = find_program('python3', required : false)
if python3.found()
    python3_config = find_program('python3-config', required : false)
    if python3_config.found()
        subdir('python')
    else
        message('python3-config not found, skipping Python wrappers.')
    endif

    java_home = run_command(python3, '-c', 'import os; print(os.environ["JAVA_HOME"])').stdout().strip()
    if java_home != ''
        subdir('java')
    else
        message('JAVA_HOME environment variable not set, skipping Java wrappers.')
    endif
else
    message('python3 not found, skipping Python and Java wrappers.')
endif
