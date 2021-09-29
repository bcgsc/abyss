[Bioinformatics Technology Lab](http://www.birollab.ca/) common code library in C++ with Python and Java wrappers.

[![Build Status](https://dev.azure.com/bcgsc/btl_public/_apis/build/status/bcgsc.btllib)](https://dev.azure.com/bcgsc/btl_public/_build/latest?definitionId=1)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/bcgsc/btllib.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/bcgsc/btllib/context:cpp)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/bcgsc/btllib.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/bcgsc/btllib/alerts/)

Platforms
---
- Linux
- MacOS

Documentation
---
[Docs page](https://bcgsc.github.io/btllib/)

Download
---
The recommended way is to download the [latest release](https://github.com/bcgsc/btllib/releases/latest).

C++
---
- Dependencies
  * GCC 5+ or Clang 4+ with OpenMP
- Copy the root `btllib` directory into your project
- Use any header from the `btllib/include` directory (pass `-I btllib/include` flag to the compiler)
- `btllib` uses `C++11` features, so that standard should be enabled at a minimum.

Python and Java
---
- Dependencies
  * GCC 5+ or Clang 4+ with OpenMP
  * Python 3.5+
  * Meson and Ninja Python3 packages, and CMake (optional -- if they are missing, they will be automatically downloaded to a temporary directory)
- Copy the root `btllib` directory into your project
- Run `btllib/compile-wrappers`
- The wrappers correspond one-to-one with C++ code so any functions and classes can be used under the same name. The only exception are nested classes which are prefixed with outer class name (e.g. `btllib::SeqReader::Flag` in C++ versus `btllib.SeqReaderFlag` in Python).
- Python
  * Use Python's `sys.path.append()` to include `btllib/python` directory
  * Include the library with `import btllib`
- Java
  * Add `btllib/java` to classpath
  * Include classes with `import btllib.*`

Contributing
---
If you want to make changes to the btllib code, first create a build directory by running `meson build` in `btllib` directory. `cd` to `build` directory and run `ninja build-sdsl` once to build the `sdsl` dependency. After that, every time you want to build the tests and wrappers, run `ninja` in `build` directory. To run the tests, run `ninja test`.

The following are the available `ninja` commands which can be run within `build` directory:
- `ninja build-sdsl` builds the sdsl-lite dependency library
- `ninja format` formats the whitespace in code (requires clang-format 8+)
- `ninja wrap` wraps C++ code for Python and Java (requires SWIG 4.0+)
- `ninja tidycheck` runs clang-tidy on C++ code and makes sure it passes (requires clang-tidy 8+)
- `ninja cppcheck` runs cppcheck on C++ code and makes sure it passes (requires cppcheck)
- `ninja` builds the tests and wrapper libraries / makes sure they compile
- `ninja test` runs the tests
- `ninja docs` generates code documentation from comments (requires Doxygen)
- `ninja complete` runs all of the above commands in the listed order

Before making a pull request, make sure to run `ninja complete` to make sure the code passes the CI test.

Credits
---
- Author: [Vladimir Nikolic](https://github.com/vlad0x00)
- Components:
  - [Hamid Mohamadi](https://github.com/mohamadi) for [ntHash](https://github.com/bcgsc/ntHash)
  - [Justin Chu](https://github.com/JustinChu) for [MIBloomFilter](https://github.com/bcgsc/btl_bloomfilter)
  - [Chase Geigle](https://github.com/skystrife) for [cpptoml](https://github.com/skystrife/cpptoml)
  - Simon Gog, Timo Beller, Alistair Moffat, and Matthias Petri for [sdsl-lite](https://github.com/simongog/sdsl-lite)