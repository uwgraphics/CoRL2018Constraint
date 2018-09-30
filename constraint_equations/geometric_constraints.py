from derive_constraint_equations import *

"""Library of constraint types. The constraint equations Phi_mat and the non-kinematic parametrization constraints 
non_kin(such as quaternion normalization) are implemented for each constraint type. """


class prismatic_constraint(constraint):


    sx, sy, sz = sym.symbols("sx sy sz")
    tx, ty, tz = sym.symbols("tx ty tz")
    wx_line, wy_line, dx, dy = sym.symbols('wx_line,wy_line,dx,dy')
    model_parameters = (sx, sy, sz, tx, ty, tz, wx_line, wy_line, dx, dy)

    def __init__(self):
        s = sym.Matrix([self.sx, self.sy, self.sz])
        t = sym.Matrix([self.tx, self.ty, self.tz])
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])
        Phi1 = ((self.r)
                - exp_map(w) * sym.Matrix([self.dx, self.dy, 0])).T * exp_map(w) * sym.Matrix([1, 0, 0])
        Phi2 = ((self.r)
                - exp_map(w) * sym.Matrix([self.dx, self.dy, 0])).T * exp_map(w) * sym.Matrix([0, 1, 0])

        Phi3 = (exp_map(w) * sym.Matrix([1, 0, 0])).T * qtoA(self.q) * s
        Phi4 = (exp_map(w) * sym.Matrix([0, 1, 0])).T * qtoA(self.q) * s
        Phi5 = (exp_map(w) * sym.Matrix([1, 0, 0])).T * qtoA(self.q) * t
        PhiS = s.T * s - sym.Matrix([1])
        PhiT = t.T * t - sym.Matrix([1])
        PhiST = s.T * t
        self.Phi_mat = Phi1.col_join(Phi2).col_join(Phi3).col_join(Phi4).col_join(Phi5)
        self.non_kin = PhiS.col_join(PhiT).col_join(PhiST)
        self.build_eqns()

class point_on_line_constraint(constraint):

    sx, sy, sz = sym.symbols("sx sy sz")
    wx_line, wy_line, dx, dy = sym.symbols('wx_line,wy_line,dx,dy')
    model_parameters = (sx, sy, sz, wx_line, wy_line, dx, dy)

    def __init__(self):
        s = sym.Matrix([self.sx, self.sy, self.sz])
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])
        Phi1 = ((self.r + qtoA(self.q) * s)
                - exp_map(w) * sym.Matrix([self.dx, self.dy, 0])).T * exp_map(w) * sym.Matrix([1, 0, 0])
        Phi2 = ((self.r + qtoA(self.q) * s)
                - exp_map(w) * sym.Matrix([self.dx, self.dy, 0])).T * exp_map(w) * sym.Matrix([0, 1, 0])

        self.Phi_mat = Phi1.col_join(Phi2)
        self.build_eqns()


class origin_on_line_constraint(constraint):

    wx_line, wy_line, dx, dy = sym.symbols('wx_line,wy_line,dx,dy')
    model_parameters = (wx_line, wy_line, dx, dy)

    def __init__(self):
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])
        Phi1 = ((self.r)
                - exp_map(w) * sym.Matrix([self.dx, self.dy, 0])).T * exp_map(w) * sym.Matrix([1, 0, 0])
        Phi2 = ((self.r)
                - exp_map(w) * sym.Matrix([self.dx, self.dy, 0])).T * exp_map(w) * sym.Matrix([0, 1, 0])

        self.Phi_mat = Phi1.col_join(Phi2)
        self.build_eqns()

class origin_on_plane_constraint(constraint):
    # (a,b,c) is the point
    # (sx,sy,sz) is the perpendicular vector on rigid body

    wx_line, wy_line, d = sym.symbols('wx_line,wy_line,d')
    model_parameters = (wx_line, wy_line, d)

    def __init__(self):
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])
        Phi1 = ((self.r)
                - exp_map(w) * sym.Matrix([0, 0, self.d])).T * exp_map(w) * sym.Matrix([0, 0, 1])
        self.Phi_mat = Phi1
        self.build_eqns()

class origin_on_arc_constraint(constraint):
    # (a,b,c) is the point
    # (sx,sy,sz) is the perpendicular vector on rigid body

    wx_line, wy_line, d,cx,cy,cz,R = sym.symbols('wx_line,wy_line,d,cx,cy,cz,R')
    model_parameters = (wx_line, wy_line, d,cx,cy,cz,R)

    def __init__(self):
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])
        c = sym.Matrix([self.cx,self.cy,self.cz])
        Phi1 = ((c)
                - exp_map(w) * sym.Matrix([0, 0, self.d])).T * exp_map(w) * sym.Matrix([0, 0, 1])
        Phi2 = ((self.r)
                - exp_map(w) * sym.Matrix([0, 0, self.d])).T * exp_map(w) * sym.Matrix([0, 0, 1])
        Phi3 = (c - self.r).T*(c - self.r)/(self.R)**2 - sym.Matrix([1.])
        Phi4 = self.R**2*0.01 # to normalize the size of R so it always finds the smallest R

        self.Phi_mat = Phi1.col_join(Phi2).col_join(Phi3)
        self.build_eqns()

class point_contact_constraint(constraint):
    # (a,b,c) is the point
    # (sx,sy,sz) is the contact location on rigid body

    a,b,c,sx,sy,sz = sym.symbols("a b c sx sy sz")
    model_parameters = (a,b,c,sx,sy,sz)

    def __init__(self):
        s = sym.Matrix([self.sx,self.sy,self.sz])
        Phi1 = sym.Matrix([self.a, self.b, self.c]) - (self.r + qtoA(self.q) * s)
        self.Phi_mat = Phi1
        self.build_eqns()


class point_on_plane_constraint(constraint):
    # (a,b,c) is the point
    # (sx,sy,sz) is the contact location on rigid body

    sx, sy, sz = sym.symbols("sx sy sz")
    wx_line, wy_line, d = sym.symbols('wx_line,wy_line,d')
    model_parameters = (sx, sy, sz, wx_line, wy_line, d)

    def __init__(self):
        s = sym.Matrix([self.sx, self.sy, self.sz])
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])
        Phi1 = ((self.r + qtoA(self.q) * s) \
                - exp_map(w) * sym.Matrix([0, 0, self.d])).T * exp_map(w) * sym.Matrix([0, 0, 1])
        self.Phi_mat = Phi1
        self.build_eqns()


class planar_constraint(constraint):
    # (a,b,c) is the point
    # (sx,sy,sz) is the perpendicular vector on rigid body

    sx, sy, sz = sym.symbols("sx sy sz")
    wx_line, wy_line, d = sym.symbols('wx_line,wy_line,d')
    model_parameters = (sx, sy, sz, wx_line, wy_line, d)

    def __init__(self):
        s = sym.Matrix([self.sx, self.sy, self.sz])
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])
        Phi1 = ((self.r)
                - exp_map(w) * sym.Matrix([0, 0, self.d])).T * exp_map(w) * sym.Matrix([0, 0, 1])
        Phi2 = (exp_map(w) * sym.Matrix([1, 0, 0])).T * qtoA(self.q) * s
        Phi3 = (exp_map(w) * sym.Matrix([0, 1, 0])).T * qtoA(self.q) * s
        PhiS = s.T * s - sym.Matrix([1])
        self.Phi_mat = Phi1.col_join(Phi2).col_join(Phi3)
        self.non_kin = PhiS
        self.build_eqns()


class axial_rotation_constraint(constraint):
    # (a,b,c) is the point
    # (sx,sy,sz) is the perpendicular vector on rigid body

    sx, sy, sz = sym.symbols("sx sy sz")
    tx, ty, tz = sym.symbols("tx ty tz")

    wx_line, wy_line, dx, dy, dz = sym.symbols('wx_line,wy_line,dx,dy,dz')
    model_parameters = (sx, sy, sz, tx, ty, tz, wx_line, wy_line, dx, dy, dz)

    def __init__(self):
        s = sym.Matrix([self.sx, self.sy, self.sz])
        t = sym.Matrix([self.tx, self.ty, self.tz])
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])

        Phi1 = self.r + qtoA(self.q) * s - exp_map(w) * sym.Matrix([self.dx, self.dy, self.dz])
        Phi2 = (exp_map(w) * sym.Matrix([0, 0, 1])).T * qtoA(self.q) * s
        Phi3 = (exp_map(w) * sym.Matrix([0, 0, 1])).T * qtoA(self.q) * t
        PhiT = t.T * t - sym.Matrix([1])
        PhiST = t.T * s

        self.Phi_mat = Phi1.col_join(Phi2).col_join(Phi3)
        self.non_kin = PhiT.col_join(PhiST)
        self.build_eqns()


class concentric_cylinder_constraint(constraint):
    # (a,b,c) is the point
    # (sx,sy,sz) is the perpendicular vector on rigid body

    sx, sy, sz = sym.symbols("sx sy sz")
    tx, ty, tz = sym.symbols("tx ty tz")

    wx_line, wy_line, dx, dy = sym.symbols('wx_line,wy_line,dx,dy')
    model_parameters = (sx, sy, sz, tx, ty, tz, wx_line, wy_line, dx, dy)

    def __init__(self):
        s = sym.Matrix([self.sx, self.sy, self.sz])
        t = sym.Matrix([self.tx, self.ty, self.tz])
        w = sym.Matrix([self.wx_line, self.wy_line, 0.0])
        regularize = 100
        Phi1 = regularize*(self.r + qtoA(self.q) * s - exp_map(w) * sym.Matrix([self.dx, self.dy, 0])).T*(exp_map(w)*sym.Matrix([1, 0, 0]))
        Phi2 = regularize*(self.r + qtoA(self.q) * s - exp_map(w) * sym.Matrix([self.dx, self.dy, 0])).T*(exp_map(w)*sym.Matrix([0, 1, 0]))

        Phi3 = (exp_map(w) * sym.Matrix([0, 0, 1])).T * qtoA(self.q) * s
        Phi4 = (exp_map(w) * sym.Matrix([0, 0, 1])).T * qtoA(self.q) * t
        PhiT = t.T * t - sym.Matrix([1])
        PhiST = t.T * s

        self.Phi_mat = Phi1.col_join(Phi2).col_join(Phi3).col_join(Phi4)
        self.non_kin = PhiT.col_join(PhiST)
        self.build_eqns()

# class orientation_constraint(constraint):
#     qsx, qsy, qsz, qsw = sym.symbols("qsx qsy qsz qsw")
#     model_parameters = (qsx, qsy, qsz, qsw)
#
#     def __init__(self):
#         qs = sym.Matrix([self.qsw,self.qsx,self.qsy,self.qsz])
#
#         Phi1 =
#
#         PhiQ = qs.T*qs - sym.Matrix([1])
#         self.Phi_mat = Phi1
#         self.non_kin = PhiQ
#         self.build_eqns()


if __name__ == "__main__":
    # c = point_contact_constraint()
    # c.build_module('point_contact_constraint','point_contact_constraint')
    #
    # c = prismatic_constraint()
    # c.build_module('prismatic_constraint','prismatic_constraint')
    #
    # c = point_on_plane_constraint()
    # c.build_module('point_on_plane_constraint','point_on_plane_constraint')
    #
    # c = planar_constraint()
    # c.build_module('planar_constraint','planar_constraint')
    #
    # c = axial_rotation_constraint()
    # c.build_module('axial_rotation_constraint','axial_rotation_constraint')
    #
    # c = concentric_cylinder_constraint()
    # c.build_module('concentric_cylinder_constraint','concentric_cylinder_constraint')
    #
    # c = point_on_line_constraint()
    # c.build_module('point_on_line_constraint','point_on_line_constraint')

    # c = origin_on_line_constraint()
    # c.build_module('origin_on_line_constraint','origin_on_line_constraint')
    #
    # c = origin_on_plane_constraint()
    # c.build_module('origin_on_plane_constraint','origin_on_plane_constraint')

    c = origin_on_arc_constraint()
    c.build_module('origin_on_arc_constraint','origin_on_arc_constraint')


