import l2g
l2g.enableDebugging()
l2g.addStreamHandler()
import l2g.equil

equilibrium = l2g.equil.getEquilibriumFromEQFile("../shot_135011_run_7_t_399.927598s.eqdsk")

from l2g.external.equilibrium_analysis import EQA

eq = EQA(equilibrium)

eq.evaluate()
print(eq.getType())
print(eq.getContactPoint())
print(eq.getLowerXPoint())
print(eq.getUpperXPoint())