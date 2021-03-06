# TrackToy
Toy MC for studying tracking resolution in precision experiments

## Installation

To build, you must have ROOT (see https://root.cern.ch/) installed, and the `bin/` directory should be on your `PATH`.
You need to have also installed and built KinKal, and source KinKal setup.sh see https://github.com/KFTrack/KinKal


There are 2 build configurations: *Debug* or *Release*.  To build

1. First, clone this repo

```bash
git clone https://github.com/KFTrack/TrackToy.git
```

2. Set up a new build directory; lets say `tracktoy_profile` for a profile build (or `tracktoy_debug` for a debug build)
```bash
mkdir tracktoy_profile
cd tracktoy_profile

```

### CMake

3. Run `cmake`, and build with `make`

```bash

cmake ../TrackToy  -DCMAKE_BUILD_TYPE=Release [Debug]

make -j <jobs to run>
```
4. Optionally, run unit tests

```bash
make test
```

Test programs will be built in the `bin/` directory. Run them with `--help` in the `build` directory to get a list of run parameters.

### Build FAQ
### Running `clang-tidy`

Make sure `clang-tidy` is installed, and add an extra argument to the `cmake` command.

```bash
mkdir tracktoy_clangtidy
cd tracktoy_clangtidy
cmake ../TrackToy -DCMAKE_BUILD_TYPE=[Release/Debug] -DENABLE_CLANG_TIDY=ON
make -j <jobs to run>
```

The code will appear to compile as normal, although in this case `clang-tidy` is actually checking the code for problems. No warnings or errors mean that no problems were found.


