@page example-apps Example Apps

\viewOnGitlab{apps/examples}

The waLBerla example applications are designed to help new users to

- learn how to create their own applications.
- discover new features.
- understand best practices for maintenance and code flow.

Examples usually consist of:

- A `CMakeLists.txt` for defining the application target and running tests.
- A `.cpp` file with the application's code, demonstrating best practices and in-code tests for the waLBerla testsuite.
- A detailed `README.md` with essential information about the example.
- Optionally:
  - A SweepGen `.py` file for generating application-specific sweeps.
  - Simulation setup parameter `.prm` files.

> **Note:** All code is subject to the GNU General Public License v3.

---

## Content

### Basic Fluid Simulations {#example-basic-fluid-simulations}

- \subpage example-DoubleShearLayer
- \subpage example-ParallelPlates
- \subpage example-FlowAroundSphere

### Framework Tools and Features {#example-tools-and-features}

- \subpage example-MeshRefinementExample
- \subpage example-TimingPool
