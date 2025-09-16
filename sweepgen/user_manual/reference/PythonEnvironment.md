# Python Environment Management

SweepGen is designed to run in a virtual Python environment.
By default, this environment is managed automatically by the waLBerla CMake system.
Refer to [](#control_managed_env) for information on controlling the base interpreter
and installed packages in the managed environment.
For advanced use cases and for development,
[](#use_external_env) details how to set up SweepGen to use an external Python environment instead.

(control_managed_env)=
## Control the Managed Environment

## Setting the Base Interpreter

WaLBerla uses [FindPython][FindPython] to locate the base Python interpreter.
Refer to its documentation for ways to affect the discovery process.

## Adding Additional Packages

To install additional packages into the SweepGen environment,
first register them using the `sweepgen_venv_require` CMake function.
This function can be invoked multiple times to add multiple requirements.
The arguments to `sweepgen_venv_require` will be directly forwarded to `pip install`,
so you can include any options `pip install` understands to affect the installation.

Calls to `sweepgen_venv_require` will only collect the set of requirements.
To perform the installation, `sweepgen_venv_populate` must be called after all
requirements are declared.

:::{card} Example

```CMake
#   First, list requirements
sweepgen_venv_require( pycowsay )  # Require a single package
sweepgen_venv_require( -r my-requirements.txt )  # Specify a requirements file

#   Then, populate the virtual environment
sweepgen_venv_populate()
```

:::

:::{error}
It is an error for your CMake system to call
`sweepgen_venv_require` after `sweepgen_venv_populate`,
or to call `sweepgen_venv_populate` more than once.
:::

[venv]: https://docs.python.org/3/library/venv.html
[FindPython]: https://cmake.org/cmake/help/latest/module/FindPython.html

(use_external_env)=
## Use an external Python Environment

To forego SweepGen's automatic environment manager,
set `WALBERLA_SWEEPGEN_MANAGED_VENV` to `OFF`.
Also, set `Python_ROOT_DIR` to point at the Python environment that should be used by
SweepGen, e.g.:

```{code-block} bash
cmake -S . -B build -DWALBERLA_SWEEPGEN_MANAGED_VENV=OFF -DPython_ROOT_DIR=<path-to-python-env>
```

You must make sure that the `sweepgen` package and its dependencies are installed
in your environment, e.g. by running `pip install -e ./sweepgen` in the waLBerla project root.
