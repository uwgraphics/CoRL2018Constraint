import numpy as np
from scipy.optimize import minimize
from robust_lsq import robust_lsq_m_estimates,robust_lsq


class InputData:
    def __init__(self,p,q,A,v,w,f,tau):
        print "initialize with arrays of each value."
        self.p = p
        self.q = q
        self.A = A
        self.v = v
        self.w = w
        self.f_global = f
        self.tau_global = tau


class ConstraintProcessing:
    c = None # constraint module
    include_velocities = True

    def __init__(self,c,include_velocities = False):
        self.c = c
        self.include_velocities = include_velocities

    def kinematic_error_per_sample(self,knowns,unknowns):
        variables = tuple(list(knowns) + list(unknowns))
        residual = np.append(self.c.Phi_mat_c(*variables), self.c.non_kin_c(*variables))
        if self.include_velocities:
            kinematics = self.c.phidelta_c(*variables)
            residual = np.append(residual, kinematics)
        return np.sum(residual**2.)

    def kinematic_lsq_error(self, r, q, v, w, unknowns, weights = None):
        if weights is None:
            weights = np.ones(len(r))
        full_objective = 0.0
        residuals = []
        for r_ind, q_ind, v_ind, w_ind, weight in zip(r, q, v, w, weights):
            knowns = list(r_ind) + list(q_ind) + list(v_ind) + list(w_ind)
            sample_error = self.kinematic_error_per_sample(knowns, unknowns)
            full_objective += np.sum(sample_error* weight)
            residuals.append(sample_error)
        return full_objective,np.array(residuals)

    def kinematic_lsq_errorJac(self,r,q,v,w,unknowns, weights = None):
        if weights is None:
            weights = np.ones(len(r))

        objective_jac_variables = np.zeros(int(self.c.len_model_parameters_c()))
        for ii, (r_ind, q_ind, v_ind, w_ind,weight) in enumerate(zip(r, q, v, w,weights)):
            knowns = list(r_ind) + list(q_ind) + list(v_ind) + list(w_ind)
            variables = tuple(list(knowns) + list(unknowns))
            objective_jac_variables += weight*(self.c.Phi_matJ_c(*variables) + \
                                        self.c.non_kinJ_c(*variables) + \
                                        self.c.phideltaJ_c(*variables))[0]
        return objective_jac_variables


    def minimize_kinematic_lsqe(self,r,q,v,w,weights,init_cond = None,numerical = False,maxiter = None):
        if init_cond == None:
            init_cond = np.random.rand(int(self.c.len_model_parameters_c()))
        objective_func = lambda X:self.kinematic_lsq_error(r,q,v,w,X,weights)[0]
        objective_jac = lambda X:self.kinematic_lsq_errorJac(r,q,v,w,X,weights)
        options = {'gtol': 1e-08, 'eps': 1e-08, 'return_all': False,
                   'maxiter': maxiter, 'norm': np.inf,'disp': True}
        if numerical == False:
            return minimize(objective_func,init_cond,
                                                  method='BFGS',jac=objective_jac,options = options)
        else:
            return minimize(objective_func,init_cond
                                                  ,method='BFGS',options=options)



    def constraint_robust_lsq(self,id,iterations=100,type = 'm_estimates'):

        X = np.concatenate((id.p, id.q, id.v, id.w), axis=1)

        def model_fit_func(X, weights):
            return self.minimize_kinematic_lsqe(X[:, :3], X[:, 3:7], X[:, 7:10], X[:, 10:], weights).x


        model_error_func = lambda X, params, weights: self.kinematic_lsq_error(X[:, :3], X[:, 3:7], X[:, 7:10],
                                                                                      X[:, 10:], params)[1] * weights
        if type == 'm_estimates':
            probs, params, errors = robust_lsq_m_estimates(model_error_func, model_fit_func, X,
                                                           iterations=iterations)
        else:
            probs, params, errors = robust_lsq(model_error_func, model_fit_func, X,
                                                           iterations=iterations)

        return probs, params, errors


    ####################################################################################################################
    ####################################################################################################################

    def force_residuals(self, id, parameters):

        tau_select = np.array([np.dot(np.transpose(A), tau) for A, tau in zip(id.A, id.tau_global)])
        residuals = []
        for p, q, v, w, f, tau in zip(id.p, id.q, id.v, id.w, id.f_global,
                                      tau_select):
            # getting reaction forces
            sol = self.get_force_torque_residuals_(p, q, f, tau, parameters)
            k = sol.x  # lagrange multipliers
            knowns = list(p) + list(q) + list(f) + list(tau)
            variables = tuple(knowns + [k.reshape(1, int(self.c.num_k_c()))] + list(parameters))
            reaction_forces = np.transpose(self.c.feq1_c(*variables))[0]
            reaction_torques = np.transpose(self.c.taueq1_c(*variables))[0]
            # getting residual forces
            residual_forces = f - reaction_forces
            residual_torques = tau - reaction_torques
            # compensating for friction
            v_hat = v / np.linalg.norm(v)
            w_hat = w / np.linalg.norm(w)
            friction_forces = np.dot(residual_forces, v_hat) * v_hat
            friction_torques = np.dot(residual_torques, w_hat) * w_hat
            f_error = np.linalg.norm(f - reaction_forces - friction_forces)
            tau_error = np.linalg.norm(tau - reaction_torques - friction_torques)
            residuals.append([f_error, tau_error])
        return residuals


    def get_force_torque_residuals_(self,r,q,fr,taur,params):
        def get_residual(r,q,fr,taur,k,params):
            knowns = list(r) + list(q) + list(fr) + list(taur)
            variables = tuple(knowns + [k.reshape(1,int(self.c.num_k_c()))] + list(params))
            residual = np.append(self.c.feq1_c(*variables), self.c.taueq1_c(*variables))
            return np.sum(residual**2)

        init_cond = np.random.rand(1,int(self.c.num_k_c()))
        objective_func = lambda K: get_residual(r, q, fr, taur,K,params)
        return minimize(objective_func, init_cond, method='BFGS')

