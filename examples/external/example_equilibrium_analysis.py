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
eq = l2g.equil.EQ(equilibrium)
eq.evaluate()

f3:float = eq.psiLCFS
f4:float = eq.psiLCFS2
t2 = eq.type_
lx2 = eq.lowXPoint
ux2 = eq.uppXPoint
Rb2, Z2, Btotal2, Bpm2 = eq.get_midplane_info()
drsep2 = eq.calculate_drsep(f3, f4)
print("Comparison between classes")
print("Type:", t1, t2)
print("Sep flux:", f1, f3)
print("2nd Sep flux:", f2, f4)
print("Drsep: ", drsep1, drsep2)
print("Lower point", lx1, lx2)
print("Upper point", ux1, ux2)
print("Rb, Z, Btotal, Bpm")
print(Rb1, Z1, Btotal1, Bpm1)
print(Rb2, Z2, Btotal2, Bpm2)
