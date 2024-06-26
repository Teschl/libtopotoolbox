# libtopotoolbox

A C++ library for the analysis of digital elevation models.

## Build instructions

From the top level directory of the repository, generate the project
buildsystem in the `build/` directory by running

```
> cmake -B build
```

and then build the library with

```
> cmake --build build
```

By default, CMake builds a static library. To build a shared library,
call

```
> cmake -B build -DBUILD_SHARED_LIBS=ON
> cmake --build build
```

## Building and running the tests

libtopotoolbox includes a test suite that can be built alongside the
library. Turn on the `TT_BUILD_TESTS` option to build the tests as well:

```
> cmake -B build -DTT_BUILD_TESTS=ON
> cmake --build build
```

Tests can then be run with

```
> cd build
> ctest
```

## Installing the library

If the library is built in the `build/` directory, it can be installed
using

```
> cmake --install build
```

This will attempt to install the library globally, which may require
administrator privileges. To install to a local path, run

```
> cmake --install build --prefix /path/to/local/installation
```

or pass the option
`-DCMAKE_INSTALL_PREFIX=/path/to/local/installation` to the initial
`cmake -B build` command.
