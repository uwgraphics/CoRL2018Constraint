import sympy as sym
import numpy as np
from scipy.optimize import minimize
from symcymod.symcymod import create_module

""" This file contains the symbolic development of the generalized constraint model.
The constraint model equations are symbolically derived using sympy and compiled into 
c/fortran code through autowrap for fast execution. Functionality to assemble and fit 
these equations using least squares regression are also provided."""

# preliminaries
def tilde(a):
    ax = a[0]
    ay = a[1]
    az = a[2]
    return sym.Matrix([[0,-az,ay],[az,0,-ax],[-ay,ax,0]])

def qtoA(q):
    e0 = q[0]
    e = sym.Matrix(q[1:4])

    A = (sym.Matrix([e0])**2 - e.T*e)[0]*sym.eye(3) + 2*e*e.T + 2*e0*tilde(e)
    return A

def getE(q):
    e0 = q[0]
    e = sym.Matrix(q[1:4])
    return (-e).row_join(tilde(e)+e0*sym.eye(3))

def getG(q):
    e0 = q[0]
    e = sym.Matrix(q[1:4])
    return (-e).row_join(-tilde(e)+e0*sym.eye(3))

def exp_map(w):
    w = sym.Matrix(w)
    # theta = ( w[0]**2 + w[1]**2 + w[2]**2 )**0.5 + 1e-30
    theta = ( w[0]**2 + w[1]**2 + w[2]**2 )**0.5
    w = w/theta
    w_hat = sym.Matrix([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
    return sym.eye(3,3) + w_hat*sym.sin(theta) + w_hat**2*(1 - sym.cos(theta))


class constraint():
    """ This is the main base class that describes the generalized model of a constraint. Different constraint model
    child classes provide Phi_mat and non_kin symbolic equations."""
    rx,ry,rz = sym.symbols("rx ry rz")
    e0,e1,e2,e3 = sym.symbols("e0 e1 e2 e3")

    frx,fry,frz = sym.symbols("frx fry frz")
    taurx,taury,taurz = sym.symbols("taurx taury taurz")

    vx,vy,vz = sym.symbols("vx vy vz")
    wx,wy,wz = sym.symbols("wx wy wz")

    q = sym.Matrix([e0,e1,e2,e3])
    r = sym.Matrix([rx,ry,rz])

    fr = sym.Matrix([frx,fry,frz])
    taur = sym.Matrix([taurx,taury,taurz])

    v = sym.Matrix([vx,vy,vz])
    w = sym.Matrix([wx,wy,wz])

    knownsF = (rx,ry,rz,e0,e1,e2,e3,frx,fry,frz,taurx,taury,taurz)
    # knownsC = (rx,ry,rz,e0,e1,e2,e3,frx,fry,frz,taurx,taury,taurz,vx,vy,vz,wx,wy,wz)
    knownsP = (rx,ry,rz,e0,e1,e2,e3,vx,vy,vz,wx,wy,wz)

    Phi_mat = None
    non_kin = sym.Matrix([0])
    taueq1 = None
    feq1 = None

    force = None
    moment = None
    work_energy = None
    model_parameters = None
    Kinematics = None
    Statics = None

    # def __init__(self,name):

    def build_eqns(self):
        Phi_mat_r = self.Phi_mat.jacobian(self.r)
        Phi_mat_q = self.Phi_mat.jacobian(self.q)
        self.num_c = len(self.Phi_mat)
        self.k = sym.MatrixSymbol('k', self.num_c,1)
        k_proxy = sym.Matrix(self.k).expand()
        # self.feq1 = Phi_mat_r.T * self.k + self.fr
        self.feq1 = Phi_mat_r.T * k_proxy + self.fr
        self.taueq1 = 1.0/2.0*getG(self.q)*Phi_mat_q.T * k_proxy + self.taur
        # position only equations
        self.phir = Phi_mat_r*self.v
        self.phiw = (1/2.0)*Phi_mat_q*(getG(self.q).T)*self.w
        self.phidelta = self.phir + self.phiw
        self.generate_jacobians()

    def generate_jacobians(self):
        # Jacobians - assuming squares are performed for least squares
        self.feq1J = (sym.Matrix(self.feq1).T * sym.Matrix(self.feq1)).jacobian(self.model_parameters)
        self.taueq1J = (sym.Matrix(self.taueq1).T*sym.Matrix(self.taueq1)).jacobian(self.model_parameters)
        self.Phi_matJ = (self.Phi_mat.T*self.Phi_mat).jacobian(self.model_parameters)
        self.non_kinJ = (self.non_kin.T*self.non_kin).jacobian(self.model_parameters)

        self.feq1Jk = (sym.Matrix(self.feq1.T*self.feq1)).jacobian(self.k)
        self.taueq1Jk = (sym.Matrix(self.taueq1.T*self.taueq1)).jacobian(self.k)
        # position only jacobians
        self.phideltaJ = (self.phidelta.T*self.phidelta).jacobian(self.model_parameters)

    def build_module(self,module_name,folder):

        # position only
        variablesP = tuple(list(self.knownsP) + list(self.model_parameters))
        variablesF = tuple(list(self.knownsF) + [self.k] + list(self.model_parameters))

        expr_tuple_list = [('Phi_mat',self.Phi_mat,variablesP),
                           ('non_kin',self.non_kin,variablesP),
                           ('phidelta', self.phidelta, variablesP),
                           ('non_kinJ', self.non_kinJ, variablesP),
                           ('Phi_matJ', self.Phi_matJ, variablesP),
                           ('phideltaJ', self.phideltaJ, variablesP),
                           ('feq1', self.feq1, variablesF),
                           ('taueq1', self.taueq1, variablesF),
                           ('feq1J', self.feq1J, variablesF),
                           ('taueq1J', self.taueq1J, variablesF),
                           ('feq1Jk', self.feq1Jk, variablesF),
                           ('taueq1Jk', self.taueq1Jk, variablesF),
                           ('len_model_parameters',len(self.model_parameters),None),
                           ('num_k',self.num_c,None)
                           ]
        create_module(module_name,expr_tuple_list,folder)


class point_contact_constraint(constraint):
    """Example constraint model"""
    # (a,b,c) is the point
    # (sx,sy,sz) is the contact location on rigid body

    a,b,c,sx,sy,sz = sym.symbols("a b c sx sy sz")
    model_parameters = (a,b,c,sx,sy,sz)

    def __init__(self):
        s = sym.Matrix([self.sx,self.sy,self.sz])
        Phi1 = sym.Matrix([self.a, self.b, self.c]) - (self.r + qtoA(self.q) * s)
        self.Phi_mat = Phi1
        self.build_eqns()


if __name__ == "__main__":
    pc = point_contact_constraint()
    pc.build_module('point_contact_constraint','point_contact_constraint')
