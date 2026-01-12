# import l2g
# l2g.enableDebugging()
# l2g.addStreamHandler()
import l2g.equil

equilibrium = l2g.equil.getEquilibriumFromEQFile("../shot_135011_run_7_t_399.927598s.eqdsk")
# print(equilibrium.grid_r.shape)
# print(equilibrium.grid_z.shape)
# print(equilibrium.psi.shape)
from l2g.external.equilibrium_analysis import EQA

eq = EQA(equilibrium)

eq.evaluate()
t1 = eq.getType()
lx1 = eq.getLowerXPoint()
ux1 = eq.getUpperXPoint()
f1 = eq.getBoundaryFluxValue()
f2 = eq.getSecondaryXFluxValue()
drsep1=eq.distanceBetweenPsiOnMidplane(f1, f2)
Rb1, Z1, Btotal1, Bpm1 = eq.getMidplaneInfo()
print("Type:", t1)
print("Sep flux:", f1)
print("2nd Sep flux:", f2)
print("Drsep: ", drsep1)
print("Lower point", lx1)
print("Upper point", ux1)
print("Rb, Z, Btotal, Bpm")
print(Rb1, Z1, Btotal1, Bpm1)
