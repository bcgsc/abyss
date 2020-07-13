BTL common code library in C++ with Python and Java wrappers.

Platforms
---
- Linux
- MacOS

Documentation
---
Open `docs/index.html` with a browser and look up any classes/functions you need.

C++
---
- Dependencies
  * GCC 4.8.1+ or Clang 3.3.0+ with OpenMP
- Copy the btllib directory in your project
- Use any header from the `btllib/include` directory

Python and Java
---
- Dependencies
  * GCC 4.8.1+ or Clang 3.3.0+ with OpenMP
  * Python 3.5+
  * Meson and Ninja Python3 packages (optional - if they are missing, they will be automatically downloaded to a temporary directory)
- Copy the btllib directory in your project
- Run `btllib/compile`
- The wrappers correspond one-to-one with C++ code so any functions and classes can be used under the same name.
- Python
  * Use Python's `sys.path.append()` to include `btllib/python` directory
  * Include the library with `import btllib`
- Java
  * Add `btllib/java` to classpath
  * Include classes with `import btllib.*`

Contributing
---
If you want to contribute code to this repo, before making a pull request, make sure to:
- Create a build directory by running `meson build` in `btllib` directory
- Run `ninja complete` in the `build` directory to generate wrappers, docs, format the code, check for any errors, etc.

`ninja complete` does the following steps, in order:
- `ninja format` formats the whitespace in code (requires clang-format 8+)
- `ninja wrap` wraps C++ code for Python and Java (requires SWIG 4.0+)
- `ninja` builds the tests and wrapper libraries / makes sure they compile
- `ninja test` runs the tests
- `ninja tidycheck` runs clang-tidy on C++ code and makes sure it passes (requires clang-tidy 8+)
- `ninja cppcheck` runs cppcheck on C++ code and makes sure it passes (requires cppcheck)
- `ninja docs` generates code documentation from comments (requires Doxygen)

Any of these can be run individually within `build` directory.
