check_LIBRARIES = libgtest.a libgtest_main.a

# gtest source files that we don't compile directly.  They are
# #included by gtest-all.cc.
EXTRA_DIST = \
  src/gtest-death-test.cc \
  src/gtest-filepath.cc \
  src/gtest-internal-inl.h \
  src/gtest-port.cc \
  src/gtest-printers.cc \
  src/gtest-test-part.cc \
  src/gtest-typed-test.cc \
  src/gtest.cc \
  include/gtest/gtest-death-test.h \
  include/gtest/gtest-message.h \
  include/gtest/gtest-param-test.h \
  include/gtest/gtest-printers.h \
  include/gtest/gtest-spi.h \
  include/gtest/gtest-test-part.h \
  include/gtest/gtest-typed-test.h \
  include/gtest/gtest.h \
  include/gtest/gtest_pred_impl.h \
  include/gtest/gtest_prod.h

libgtest_a_CPPFLAGS = -I$(top_srcdir) -isystem include
libgtest_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-missing-field-initializers
libgtest_a_SOURCES = src/gtest-all.cc

libgtest_main_a_CPPFLAGS = -I$(top_srcdir) -isystem include
libgtest_main_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-missing-field-initializers
libgtest_main_a_SOURCES = src/gtest_main.cc src/gtest-all.cc
