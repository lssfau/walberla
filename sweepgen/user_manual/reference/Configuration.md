# Configuration

```{eval-rst}
.. module:: sweepgen.build_config
```

## Build Configuration

SweepGen learns its configuration from the waLBerla build system
and passes it on to the pystencils code generation system.
This includes information about the target architecture (CPU, CUDA, HIP, ...)
and optimization options (OpenMP, vectorization).
For all options passed this way, refer to {any}`WalberlaBuildConfig`.

The configuration object controlling a running generator script
is maintained in the `pystencils-sfg` session context.
You can get a reference to the configuration using {any}`WalberlaBuildConfig.from_sfg`
and use it to override certain configuration options,
such as the code generation {any}`target <WalberlaBuildConfig.target>`.

```{eval-rst}
.. autoclass:: WalberlaBuildConfig
    :members:

.. autoclass:: ConfigOverrides
    :members:
```

## Mock the Build System for Debugging

When run outside of waLBerla's build system for debugging,
SweepGen falls back to a mockup build configuration set through the
{any}`sweepgen.build_config.DEBUG` class:

```{eval-rst}
.. autoclass:: DEBUG
    :members:
```
