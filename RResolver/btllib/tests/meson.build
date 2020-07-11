tests_files = run_command('../scripts/get_tests_files').stdout().strip().split()

foreach file : tests_files
    t = file.split('.cpp')[0].split('.cxx')[0].split('/')[-1]
    test(t, executable(t, file, include_directories : btllib_include))
endforeach