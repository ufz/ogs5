[![Tag](https://img.shields.io/github/tag/ufz/ogs5.svg?style=flat-square)](https://github.com/ufz/ogs5/releases)
[![BSD License (modified)](http://img.shields.io/badge/license-BSD-blue.svg?style=flat-square)](https://github.com/ufz/ogs5/blob/master/LICENSE.txt)
[![Travis](https://img.shields.io/travis/ufz/ogs5.svg?style=flat-square)](https://travis-ci.org/ufz/ogs5)

# OGS-5

- General homepage: https://www.opengeosys.org/ogs-5
- Keyword documentation: https://ogs5-keywords.netlify.com
- Benchmarks: https://github.com/ufz/ogs5-benchmarks
- [ogs-users](https://groups.google.com/forum/#!forum/ogs-users) mailing list

## Quickstart option 1: binary download

- Go to the [GitHub release page](https://github.com/ufz/ogs5/releases)
- Download the appropriate archive (`*.zip` for Windows, `*.tar.gz` for Linux, `*.dmg` for macOS) and unpack it
- You will find the `ogs`-binary inside the `bin`-folder
- *Optional:* [Download](https://github.com/ufz/ogs5/releases/tag/data-explorer-5) and unpack the Data Explorer for OGS-5 too
- [Download](https://github.com/ufz/ogs5-benchmarks/archive/master.zip) and unpack benchmark files

## Quickstart option 2: build from source

### Prerequisites

- Install a compiler
  - Win: Download and install [Visual Studio Community 2017](https://www.visualstudio.com/de/thank-you-downloading-visual-studio/?sku=Community&rel=15), during installation select the *workload* `Desktop Development with C++`, uncheck everything else
  - Linux: Install the following packages: `make` and `build-essential`
  - macOS: Install Xcode from the AppStore and run the following in the command line: `xcode-select --install`
- Install [Git](https://git-scm.com)
- Install [CMake](https://cmake.org/download/)

### Source code

Get the source code either by [downloading as a zip-file](https://github.com/ufz/ogs5/archive/master.zip) or with Git:

```bash
git clone https://github.com/ufz/ogs5.git
```

### Build

Then set-up a build directory, configure your build and compile the code:

```bash
cd [source-directory]
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release # Configures your build
cmake --build . --config Release    # Compiles the code
```

### Further info

- [How to Build a CMake-Based Project](http://preshing.com/20170511/how-to-build-a-cmake-based-project/), a nice overview how CMake works, contains also information for building with [Visual Studio](http://preshing.com/20170511/how-to-build-a-cmake-based-project/#building-with-visual-studio)
- [Running CMake](https://cmake.org/runningcmake/) describes various to run CMake also from a GUI and how to set configuration options
- [Git tutorial](https://www.atlassian.com/git/tutorials) and [Git FAQs](https://github.com/k88hudson/git-flight-rules)

## Quickstart: Run an OGS benchmark

Open a command prompt and run the following:

```bash
cd [benchmark-directory] # e.g. ogs5-benchmarks-master/H/Theis/GWF_Theis_2D
[path-to-ogs-exe-folder]/ogs [benchmark name] # e.g. ../../build/bin/ogs GWF_Theis_2d
```

## Contributing: Using Git and GitHub

To implement new features, every developer
1. [forks](https://help.github.com/articles/fork-a-repo/) this repository to have their own repo
2. creates a [branch](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging) from master and implement their stuff
3. [pushes](https://help.github.com/articles/pushing-to-a-remote/) the local branch to its GitHub repository
4. makes a [pull request](https://help.github.com/articles/creating-a-pull-request/) from its GitHub repository to the `master` branch in the `ufz/ogs5` repository

