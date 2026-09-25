# -*- coding: utf-8 -*-
import os, sys, time
import openmm
print("openmm", openmm.__version__, flush=True)
from openmm import Platform, Context, System, NonbondedForce, LangevinIntegrator, unit

plat = Platform.getPlatformByName("CUDA")
print("plataforma CUDA cargada, propiedades:", plat.getPropertyValue("CudaDevice") if False else "ok", flush=True)
t0 = time.time()
s = System()
s.addParticle(12.0)
f = NonbondedForce()
f.addParticle(0.0, 0.3, 0.0)
f.setNonbondedMethod(NonbondedForce.NoCutoff)
s.addForce(f)
integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)
ctx = Context(s, integ, plat)
ctx.setPositions([[0, 0, 0]] * unit.nanometer)
e = ctx.getState(getEnergy=True).getPotentialEnergy()
print("CUDA energia OK:", e, "| init %.1fs" % (time.time() - t0), flush=True)
print("dispositivo:", plat.getPropertyValue("CudaDevice"), plat.getPropertyValue("UseCpuPme"), flush=True)
