import json
import yaml
from jsonschema import validate, ValidationError

# Additional schema checks can be added here as file, schema tuples
validate_files = [("CITATION.cff", "utilities/schema-validation/.cffschema")]

success = True

for file, schema_file in validate_files:
    print(f"Validating {file} against {schema_file}...")
    # Load schema
    with open(schema_file, "r") as f:
        schema = json.load(f)

    # Load and validate CITATION.cff
    with open(file, "r") as f:
        cff_data = yaml.safe_load(f)

    try:
        validate(instance=cff_data, schema=schema)
        print(f"✓ {file} is valid according to {schema_file}")
    except ValidationError as e:
        print(f"✗ {file} validation failed: {e.message}")
        success = False

if not success:
    exit(1)
