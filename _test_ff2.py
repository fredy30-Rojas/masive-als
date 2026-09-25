import openff.toolkit
print("toolkit path:", openff.toolkit.__file__)
try:
    from openff.toolkit.typing.engines.smirnoff.forcefield import get_available_force_fields
    print("AVAILABLE:", get_available_force_fields())
except Exception as e:
    print("get_available err:", e)
# buscar la lista interna de paths
try:
    from openff.toolkit.typing.engines.smirnoff import forcefield as ffmod
    import inspect
    src = inspect.getsource(ffmod._get_forcefield_data_paths)
    print(src[:1500])
except Exception as e:
    print("path err:", e)
import openff.toolkit.data
print("data dir:", openff.toolkit.data.__path__)