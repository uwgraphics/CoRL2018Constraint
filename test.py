from fit_constraints.fit_constraints import ConstraintProcessing, InputData
import point_contact_constraint.point_contact_constraint as pcc
import numpy as np

pccProc = ConstraintProcessing(pcc,include_velocities=False)

p = np.array([[1.0, 0, 0], [0, 1, 0]])
q = np.array([[1, 0, 0, 0], [0, 1.0 / 2 ** 0.5, 1.0 / 2 ** 0.5, 0]])
A = [np.eye(3),np.eye(3)]
f = np.array([[0., 1, 0], [1, 0, 0]])
tau = np.array([[0., 0, -1], [0, 0, -1]])
v = np.array([[0.001, 0., 0], [-0.001, 0, 0.0]])
w = np.array([[100.00, 0, 0], [0.001, 0, -0]])

id = InputData(p,q,A,v,w,f,tau)

anss = pccProc.constraint_robust_lsq(id)
print anss
test_inputs = tuple(list(np.random.rand(13)) + [np.random.rand(1,3)] + list(np.random.rand(6)))
print pccProc.force_residuals(id,anss[1])[0]

