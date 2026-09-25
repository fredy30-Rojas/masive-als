from openff.toolkit import ForceField
p = r"C:\Users\Fredy\masive-als\rescoring_env\Lib\site-packages\openff\toolkit\data\forcefield\openff-2.1.0.offxml"
ff = ForceField(p)
print("ForceField por ruta OK, handlers:", len(ff.registered_parameter_handlers))