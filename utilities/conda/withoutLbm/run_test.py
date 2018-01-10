import waLBerla

# Check that C++ exports are available
assert waLBerla.cpp_available

# Test calling of a function taking a string (fails if there is a ABI compatibility issue)
waLBerla.log_devel("Test successful")
