# FAQ: Common Problems When Building waLBerla  {#faq}

\tableofcontents

## CMake

### I am working on an HPC cluster with an outdated version of CMake. Can I still build waLBerla?  {#faq-outdated-cmake}

Yes; there are several ways to get an up-to-date version of CMake even without super-user rights on a machine.

#### Install CMake using pip (recommended)

If you have Python available, you may create and activate a new
[virtual environment](https://docs.python.org/3/library/venv.html)
and install CMake into it using `pip`:

```bash
pip install cmake
```

You can ensure that you are using the `pip`-installed CMake version like this:

```bash
which cmake  # should point to a path inside your virtual environment
cmake --version   # should be the latest CMake release
```

#### Download a Binary Distribution

You can download a binary distribution of CMake directly from [cmake.org](https://cmake.org/download/)
and install it into your user directory.
In this case, you might have to manually add the installation directory to `PATH`.
