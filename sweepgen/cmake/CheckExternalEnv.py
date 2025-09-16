import sys


def error(msg: str):
    msg = "Required Python packages could not be found:\n" + msg
    print(msg, file=sys.stderr)
    sys.exit(1)

try:
    import sweepgen
except ImportError:
    error("sweepgen is not installed in the current Python environment.")

print("Found required Python packages.", end="")
sys.exit(0)
