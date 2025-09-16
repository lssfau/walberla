# IDE Configuration & Debugging

## IDE Configuration

When writing generator scripts in an IDE,
make sure to point your IDE to the correct Python interpreter to
get the most out of syntax highlighting and linting.
By default, the SweepGen Python interpreter is located at
`<walberla-build-directory>/sweepgen-venv/bin/python`.
If you are using an [external Python environment](#use_external_env),
select its interpreter instead.

## Debugging Generator Scripts

Since SweepGen generator scripts are designed to be called by and communicate with
waLBerla's build system, debugging them is not completely straightforward, but also
not difficult.
To run and debug a script outside of the build system,
one merely needs to mimic the environment otherwise provided by waLBerla.

## Activate the Python Environment

When working in a console, activate the SweepGen Python environment:

```{code-block} bash
source <walberla-build-directory>/sweepgen-venv/bin/activate
```

## Mimic the Build Environment

To mimic a specific build environment, add the following lines at the top of your generator script:

```{code-block} python
from sweepgen.build_config import DEBUG
DEBUG.use_cpu_default()
```

These two lines instruct SweepGen to fall back to a default CPU build configuration
when it cannot get a configuration from waLBerla.
Refer to the documentation of {any}`build_config.DEBUG` for other available configurations,
e.g. using CUDA or HIP.

## Run the Generator Script

For debugging, run your generator script using the following command:

```{code-block} bash
python [your-script].py --sfg-output-dir out
```

This will emit all generated files to the `out` directory,
which you might want to add to your `.gitignore`.
