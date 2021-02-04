#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from copy import deepcopy
import numpy as np
from scipy import sparse as sp
from scipy.linalg import eig
from scipy.sparse.linalg import spsolve

class regals:

    def __init__(self, I, err):
        self.I = I
        self.err = err

    def fit_concentrations(self, mix):

        mix = deepcopy(mix)

        H = mix.H_concentration
        [AA, Ab] = mix.concentration_problem(self.I, self.err)

        u = spsolve(AA + H, Ab)

        u = np.split(u, np.cumsum(mix.k_concentration)[:-1])

        mix.u_concentration = u
        return mix

    def fit_profiles(self, mix):

        mix = deepcopy(mix)

        H = mix.H_profile
        [AA, Ab] = mix.profile_problem(self.I, self.err)

        u = spsolve(AA + H, Ab)

        u = np.split(u, np.cumsum(mix.k_profile)[:-1])

        n = mix.norm_profile
        u = [uj / nj if nj != 0 else uj for uj, nj in zip(u, n)]

        mix.u_profile = u
        return mix

    def step(self, mix):

        new_mix = self.fit_concentrations(self.fit_profiles(mix));

        resid = (self.I - new_mix.I_reg) / self.err

        params = {}
        params['x2'] = np.mean(resid ** 2)
        params['delta_concentration'] = np.sum(np.abs(new_mix.concentrations - mix.concentrations),0)
        params['delta_profile'] = np.sum(np.abs(new_mix.profiles - mix.profiles),0)
        params['delta_u_concentration'] = np.array([np.sum(np.abs(nupk - upk)) for nupk, upk in zip(new_mix.u_concentration, mix.u_concentration)])
        params['delta_u_profile'] = np.array([np.sum(np.abs(nupr - upr)) for nupr, upr in zip(new_mix.u_profile, mix.u_profile)])

        return [new_mix, params, resid]

    def run(self, mix, stop_fun = None, update_fun = None):

        if stop_fun is None:
            stop_fun = lambda num_iter, params: [num_iter >= 10, 'max_iter']

        if update_fun is None:
            update_fun = lambda num_iter, new_mix, params, resid: True

        num_iter = 0
        while True:
            num_iter += 1
            [mix, params, resid] = self.step(mix)

            update_fun(num_iter, mix, params, resid)

            [if_exit, exit_cond] = stop_fun(num_iter, params)
            if if_exit:
                break

        return [mix, params, resid, exit_cond]



class mixture:

    def __init__(self, components, lambda_concentration = np.array([]), lambda_profile = np.array([]), u_concentration = [], u_profile = []):
        self.components = components
        self._lambda_concentration = lambda_concentration
        self._lambda_profile = lambda_profile
        self.u_concentration = u_concentration
        self.u_profile = u_profile

        self.Nc = len(self.components)

        self.Nx = self.components[0].concentration.Nx
        self.Nq = self.components[0].profile.Nq
        self.k_concentration = np.array([comp.concentration.k for comp in self.components])
        self.k_profile = np.array([comp.profile.k for comp in self.components])

        if len(self.u_concentration) == 0:
            self.u_concentration = [comp.concentration.u0 for comp in components]

        if len(self.u_profile) == 0:
            self.u_profile = [comp.profile.u0 for comp in components]

        if len(self._lambda_concentration) == 0:
            self._lambda_concentration = np.zeros(self.Nc)

        if len(self._lambda_profile) == 0:
            self._lambda_profile = np.zeros(self.Nc)

        self._calc_params()

    @property
    def lambda_concentration(self):
        return self._lambda_concentration

    @lambda_concentration.setter
    def lambda_concentration(self, val):
        if len(val) == 0:
            self._lambda_concentration = np.zeros(self.Nc)
        else:
            self._lambda_concentration = val
            self._calc_params()

    @property
    def lambda_profile(self):
        return self._lambda_profile

    @lambda_profile.setter
    def lambda_profile(self, val):
        if len(val) == 0:
            self._lambda_profile = np.zeros(self.Nc)
        else:
            self._lambda_profile = val
            self._calc_params()

    @property
    def I_reg(self):
        return self.profiles @ self.concentrations.T

    def estimate_concentration_lambda(self,err,ng):
        AA = self.concentration_problem(np.zeros((self.Nq,self.Nx)),err,False).todense()

        split_pos = np.cumsum(self.k_concentration)[:-1]
        AA = [np.hsplit(AAi, split_pos) for AAi in np.vsplit(AA, split_pos)]

        ll = np.zeros(self.Nc)
        L = [comp.concentration.L for comp in self.components]
        for k in range(self.Nc):
            d = np.real(eig(AA[k][k], (L[k].T @ L[k]).todense())[0])
            ll[k] = mixture.ng2lambda(d, ng[k])

        return ll

    def estimate_profile_lambda(self,err,ng):
        AA = self.profile_problem(np.zeros((self.Nq,self.Nx)),err,False).todense()

        split_pos = np.cumsum(self.k_profile)[:-1]
        AA = [np.hsplit(AAi, split_pos) for AAi in np.vsplit(AA, split_pos)]

        ll = np.zeros(self.Nc)
        L = [comp.profile.L for comp in self.components]
        for k in range(self.Nc):
            d = np.real(eig(AA[k][k], (L[k].T @ L[k]).todense())[0])
            ll[k] = mixture.ng2lambda(d, ng[k])

        return ll

    def concentration_problem(self, I, err, calc_Ab = True):

        w = 1 / np.mean(err,1)

        A = [comp.concentration.A for comp in self.components]

        D = w[:,np.newaxis] * I
        y = self.profiles
        y = w[:,np.newaxis] * y

        AA = [[(y[:,k1] @ y[:,k2]) * (A[k1].T @ A[k2]) for k2 in range(self.Nc)] for k1 in range(self.Nc)]
        AA = sp.vstack((sp.hstack(tuple(AAi)) for AAi in AA))

        if calc_Ab == True:
            Ab = [A[k].T @ (D.T @ y[:,k]) for k in range(self.Nc)]
            Ab = np.hstack(tuple(Ab))
            return [AA, Ab]
        else:
            return AA

    def profile_problem(self, I, err, calc_Ab = True):

        w = 1 / np.mean(err,1)

        A = [comp.profile.A for comp in self.components]
        A = [sp.diags(w,0) @ Ai for Ai in A]

        D = w[:,np.newaxis] * I
        c = self.concentrations

        AA = [[sp.csr_matrix((c[:,k1] @ c[:,k2]) * (A[k1].T @ A[k2])) for k2 in range(self.Nc)] for k1 in range(self.Nc)]
        AA = sp.vstack((sp.hstack(tuple(AAi)) for AAi in AA))

        if calc_Ab == True:
            Ab = [A[k].T @ (D @ c[:,k]) for k in range(self.Nc)]
            Ab = np.hstack(tuple(Ab))
            return [AA, Ab]
        else:
            return AA

    def extract_concentration(self,I,err,k):

        notk = np.setdiff1d(np.arange(self.Nc), k)
        c = self.concentrations
        y = self.profiles
        D = I - y[:,notk] @ c[:,notk].T
        yk = y[:,k]
        w = 1/np.mean(err,1)
        m = (w**2 * yk) / (yk @ (w**2 * yk))
        pk = D.T @ m
        sigmak = np.sqrt(err.T**2 @ m**2)

        return [pk, sigmak]

    def extract_profile(self,I,err,k):

        notk = np.setdiff1d(np.arange(self.Nc), k)
        c = self.concentrations
        y = self.profiles
        D = I - y[:,notk] @ c[:,notk].T
        ck = c[:,k]
        m = ck/(ck @ ck)
        Ik = D @ m
        sigmak = np.sqrt(err**2 @ m**2)

        return [Ik, sigmak]

    @property
    def concentrations(self):
        return np.hstack(tuple(comp.concentration.A @ upk[:,np.newaxis] for comp, upk in zip(self.components, self.u_concentration)))

    @property
    def profiles(self):
        return np.hstack(tuple(comp.profile.A @ upr[:,np.newaxis] for comp, upr in zip(self.components, self.u_profile)))

    @property
    def norm_concentration(self):
        return [comp.concentration.norm(upk) for comp, upk in zip(self.components, self.u_concentration)]

    @property
    def norm_profile(self):
        return [comp.profile.norm(upr) for comp, upr in zip(self.components, self.u_profile)]

    @staticmethod
    def ng2lambda(dd, ng):

        ng0 = np.count_nonzero(np.isinf(dd))

        dd = dd[np.logical_and(dd >= 0,~np.isinf(dd))]

        lambda_list = np.logspace(np.amax(np.log10(dd)) + 2, np.amin(np.log10(dd)) - 2, 51)

        ng_list = np.zeros(len(lambda_list))
        for j in range(len(ng_list)):
            ng_list[j] = ng0 + sum(dd / (dd + lambda_list[j]))

        if ng < ng_list[0]:
            optimal_lambda = np.inf
        elif ng > ng_list[-1]:
            optimal_lambda = 0
        else:
            optimal_lambda = 10 ** np.interp(ng, ng_list, np.log10(lambda_list))

        return optimal_lambda

    def _calc_params(self):


        conc_L = [comp.concentration.L for comp in self.components]
        conc_B = [Lpk * (lbdpk ** 0.5) for Lpk, lbdpk in zip(conc_L, self.lambda_concentration)]
        conc_B = sp.block_diag(conc_B)
        self.H_concentration = conc_B.T @ conc_B


        prof_L = [comp.profile.L for comp in self.components]
        prof_B = [Lpr * (lbdpr ** 0.5) for Lpr, lbdpr in zip(prof_L, self.lambda_profile)]
        prof_B = sp.block_diag(prof_B)
        self.H_profile = prof_B.T @ prof_B

class component:

    def __init__(self, concentration, profile):
        self.concentration = concentration
        self.profile = profile



class concentration_class:

    def __init__(self, concentration_type, *arg, **kwarg):
        self.concentration_type = concentration_type

        if self.concentration_type == 'simple':
            self._regularizer = concentration_simple(*arg, **kwarg)

        elif self.concentration_type == 'smooth':
            self._regularizer = concentration_smooth(*arg, **kwarg)

        else:
            raise ValueError('unexpected concentration type')

    def __getattr__(self,attr):
        return super().__getattribute__('_regularizer').__getattribute__(attr) #super() usage for deepcopy



class concentration_simple:

    def __init__(self, x, xmin, xmax):
        self.x = x
        self._xmin = xmin
        self._xmax = xmax

        self._calc_params()

    @property
    def xmin(self):
        return self._xmin

    @xmin.setter
    def xmin(self, val):
        self._xmin = val
        self._calc_params()

    @property
    def xmax(self):
        return self._xmax

    @xmax.setter
    def xmax(self, val):
        self._xmax = val
        self._calc_params()

    def norm(self,u):
        return np.mean(u**2)**0.5

    def _calc_params(self):
        self.Nx = len(self.x)

        # Calculate F
        is_in_concentration = np.isin(self.x,self.w)
        ind_in_w = np.nonzero(self.x[:,np.newaxis] == self.w)[1] #assumes no repetition in x
        v = np.arange(self.Nx)
        self.F = sp.csr_matrix((np.ones(self.Nw), (v[is_in_concentration], ind_in_w)), shape=(self.Nx, self.Nw))

        self.w = self.x[np.logical_and(self.x >= self.xmin, self.x <= self.xmax)] #assumes no repetition in x
        self.Nw = len(self.w)
        self.k = self.Nw
        self.u0 = np.ones(self.k)
        self.L = sp.eye(self.Nw)
        self.A = self.F
        self.y0 = self.A @ self.u0



class concentration_smooth:

    def __init__(self, x, xmin, xmax, Nw = 50, is_zero_at_xmin = True, is_zero_at_xmax = True):
        #Changed input parameter order because python requires default parameters to follow non-default ones
        self.x = x
        self._xmin = xmin
        self._xmax = xmax
        self._Nw = Nw
        self._is_zero_at_xmin = is_zero_at_xmin
        self._is_zero_at_xmax = is_zero_at_xmax

        self._calc_params()

    @property
    def xmin(self):
        return self._xmin

    @xmin.setter
    def xmin(self, val):
        self._xmin = val
        self._calc_params()

    @property
    def xmax(self):
        return self._xmax

    @xmax.setter
    def xmax(self, val):
        self._xmax = val
        self._calc_params()

    @property
    def Nw(self):
        return self._Nw

    @Nw.setter
    def Nw(self, val):
        self._Nw = val
        self._calc_params()

    @property
    def is_zero_at_xmin(self):
        return self._is_zero_at_xmin

    @is_zero_at_xmin.setter
    def is_zero_at_xmin(self, val):
        self._is_zero_at_xmin = val
        self._calc_params()

    @property
    def is_zero_at_xmax(self):
        return self._is_zero_at_xmax

    @is_zero_at_xmax.setter
    def is_zero_at_xmax(self, val):
        self._is_zero_at_xmax = val
        self._calc_params()

    def norm(self,u):
        return np.mean(u**2)**0.5

    def _calc_params(self):
        self.Nx = len(self.x)
        self.k = self.Nw - self.is_zero_at_xmin - self.is_zero_at_xmax
        self.dw = (self.xmax - self.xmin)/(self.Nw-1)
        self.w = np.linspace(self.xmin, self.xmax, self.Nw)

        # Calculate F
        ix = np.searchsorted(self.w,self.x,side='right')
        ix = ix - 1
        is_in_concentration = np.logical_and(ix > -1, ix < self.Nw - 1)
        v = np.arange(self.Nx)
        ix = ix[is_in_concentration]
        v = v[is_in_concentration]
        u = (self.x[is_in_concentration] - self.w[ix]) * (1/self.dw)
        self.F = sp.csr_matrix((np.concatenate((1-u,u)), (np.concatenate((v,v)), np.concatenate((ix,ix+1)))),shape = (self.Nx, self.Nw))

        # Calculate u0
        if self.is_zero_at_xmax and self.is_zero_at_xmin:
            u0_tmp =  1 - ( 2*(self.w-self.xmin)/(self.xmax-self.xmin) - 1) ** 2
            self.u0 = u0_tmp[1:-1]
        elif self.is_zero_at_xmax and (not self.is_zero_at_xmin):
            u0_tmp = (self.xmax-self.w)/(self.xmax-self.xmin)
            self.u0 = u0_tmp[:-1]
        elif (not self.is_zero_at_xmax) and self.is_zero_at_xmin:
            u0_tmp = (self.w-self.xmin)/(self.xmax-self.xmin)
            self.u0 = u0_tmp[1:]
        else:
            self.u0 = np.ones(self.Nw)

        # Calculate L
        L_tmp = sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(0,self.Nw-2))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(1,self.Nw-1))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(2,self.Nw))), shape=(self.Nw-2, self.Nw))
        if self.is_zero_at_xmax:
            L_tmp = L_tmp[:,:-1]
        if self.is_zero_at_xmin:
            L_tmp = L_tmp[:,1:]
        self.L = L_tmp

        #Calculate A
        A_tmp = self.F
        if self.is_zero_at_xmax:
            A_tmp = A_tmp[:,:-1]
        if self.is_zero_at_xmin:
            A_tmp = A_tmp[:,1:]
        self.A = A_tmp

        self.y0 = self.A @ self.u0



class profile_class:

    def __init__(self, profile_type, *arg, **kwarg):
        self.profile_type = profile_type

        if self.profile_type == 'simple':
            self._regularizer = profile_simple(*arg, **kwarg)

        elif self.profile_type == 'smooth':
            self._regularizer = profile_smooth(*arg, **kwarg)

        elif self.profile_type == 'realspace':
            self._regularizer = profile_real_space(*arg, **kwarg)

        else:
            raise ValueError('unexpected profile type')

    def __getattr__(self,attr):
        return super().__getattribute__('_regularizer').__getattribute__(attr) #super() usage for deepcopy



class profile_simple:

    def __init__(self, q):
        self.q = q

        self.Nq = len(self.q)
        self.F = sp.eye(self.Nq)
        self.k = self.Nq
        self.u0 = np.ones(self.k)
        self.w = np.ones(self.k)
        self.Nw = self.Nq
        self.L = sp.eye(self.Nq)
        self.A = self.F
        self.y0 = self.A @ self.u0

    def norm(self,u):
        return np.mean(u**2)**0.5



class profile_smooth:

    def __init__(self, q, Nw = 50):
        self.q = q
        self._Nw = Nw

        self._calc_params()

    @property
    def Nw(self):
        return self._Nw

    @Nw.setter
    def Nw(self, val):
        self._Nw = val
        self._calc_params()

    def norm(self,u):
        return np.mean(u**2)**0.5

    def _calc_params(self):
        self.Nq = len(self.q)
        self.k = self.Nw
        self.u0 = np.ones(self.k)
        self.qmin = np.amin(self.q)
        self.qmax = np.amax(self.q)
        self.dw = (self.qmax - self.qmin) / (self.Nw - 1)
        self.w = np.linspace(self.qmin, self.qmax, self.Nw)

        # Calculate F
        ix = np.searchsorted(self.w,self.q,side='right')
        ix = ix - 1
        ix[ix == -1] = 0
        ix[ix == self.Nw - 1] = self.Nw - 2
        v = np.arange(self.Nq)
        u = (self.q - self.w[ix]) * (1/self.dw)
        self.F = sp.csr_matrix((np.concatenate((1-u,u)), (np.concatenate((v,v)), np.concatenate((ix,ix+1)))),shape = (self.Nq, self.Nw))

        self.A = self.F

        # Calculate L
        self.L = sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(0,self.Nw-2))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(1,self.Nw-1))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(2,self.Nw))), shape=(self.Nw-2, self.Nw))



class profile_real_space:

    def __init__(self, q, dmax, Nw = 50, is_zero_at_r0 = True, is_zero_at_dmax = True):
        self.q = q
        self._dmax = dmax
        self._Nw = Nw
        self._is_zero_at_r0 = is_zero_at_r0
        self._is_zero_at_dmax = is_zero_at_dmax

        self._calc_params()

    @property
    def dmax(self):
        return self._dmax

    @dmax.setter
    def dmax(self, val):
        self._dmax = val
        self._calc_params()

    @property
    def Nw(self):
        return self._Nw

    @Nw.setter
    def Nw(self, val):
        self._Nw = val
        self._calc_params()

    @property
    def is_zero_at_r0(self):
        return self._is_zero_at_r0

    @is_zero_at_r0.setter
    def is_zero_at_r0(self, val):
        self._is_zero_at_r0 = val
        self._calc_params()

    @property
    def is_zero_at_dmax(self):
        return self._is_zero_at_dmax

    @is_zero_at_dmax.setter
    def is_zero_at_dmax(self, val):
        self._is_zero_at_dmax = val
        self._calc_params()

    def norm(self,u):
        weight = 4 * np.pi * self.dw * np.ones(self.Nw)
        weight[[0,-1]] = 0.5* weight[[0,-1]]
        if self.is_zero_at_dmax:
            weight = weight[:-1]
        if self.is_zero_at_r0:
            weight = weight[1:]
        return np.sum(weight * u)

    def _calc_params(self):
        self.Nq = len(self.q)

        self.k = self.Nw - int(self.is_zero_at_dmax) - int(self.is_zero_at_r0)
        self.dw = self.dmax / (self.Nw - 1)
        self.w = self.dmax * np.linspace(0,1,self.Nw)

        # Calculate u0
        u0_tmp = 1 - (2 * self.w / self.dmax - 1) ** 2
        u0_tmp = u0_tmp / (4 * np.pi * self.dw * np.sum(u0_tmp))
        if self.is_zero_at_dmax:
            u0_tmp = u0_tmp[:-1]
        if self.is_zero_at_r0:
            u0_tmp = u0_tmp[1:]
        self.u0 = u0_tmp

        # Calculate F
        F_tmp = 4 * np.pi * self.dw * np.sinc(np.outer(self.q,self.w) / np.pi)
        F_tmp[:,[0,-1]] = 0.5*F_tmp[:,[0,-1]]
        self.F = F_tmp

        # Calculate A
        A_tmp = self.F
        if self.is_zero_at_dmax:
            A_tmp = A_tmp[:,:-1]
        if self.is_zero_at_r0:
            A_tmp = A_tmp[:,1:]
        self.A = A_tmp

        # Calculate L
        L_tmp = sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(0,self.Nw-2))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(1,self.Nw-1))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(2,self.Nw))), shape=(self.Nw-2, self.Nw))
        if self.is_zero_at_dmax:
            L_tmp = L_tmp[:,:-1]
        if self.is_zero_at_r0:
            L_tmp = L_tmp[:,1:]
        self.L = L_tmp

        self.y0 = self.A @ self.u0
