from openff.toolkit.typing.engines.smirnoff import get_available_force_fields
print("AVAILABLE:", get_available_force_fields())
from openff.toolkit import ForceField
ff = ForceField("openff-2.1.0")
print("ForceField OK, handlers:", len(ff.registered_parameter_handlers))