import emcee
import ctypes
import os
import numpy as np
import sys
import copy

dirname = os.path.split(__file__)[0]
lib_name = os.path.join(dirname, "./like_mix_6x2pt.so")
lib=ctypes.cdll.LoadLibrary(lib_name)
double = ctypes.c_double

Double10 = double*10

write_cosmolike_datavector = lib.write_datavector_wrapper


initcosmo=lib.init_cosmo_runmode
initcosmo.argtypes=[ctypes.c_char_p]

initbins=lib.init_binning_fourier
initbins.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double]

initscalecuts=lib.init_scalecuts
initscalecuts.argtypes=[ctypes.c_double,ctypes.c_double]

initsources=lib.init_source_sample_mpp
initsources.argtypes=[ctypes.c_char_p,ctypes.c_int]

initlenses=lib.init_lens_sample_mpp
initlenses.argtypes=[ctypes.c_char_p,ctypes.c_int,Double10, Double10,ctypes.c_double]

initia=lib.init_IA_mpp
initia.argtypes=[ctypes.c_int]

initthetas=lib.init_sample_theta_s
initthetas.argtypes=[]

initprobes=lib.init_probes
initprobes.argtypes=[ctypes.c_char_p]

initdata_fourier=lib.init_data_fourier
initdata_fourier.argtypes=[ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
#init_filenames=lib.init_filenames
#init_filenames.argtypes=[ctypes.c_char_p, ctypes.c_char_p]

setprior_m=lib.set_shear_priors_mpp
setprior_m.argtypes =[Double10,Double10]

setprior_wlphotoz=lib.set_wlphotoz_priors_mpp
setprior_wlphotoz.argtypes =[Double10,Double10]

setprior_clusteringphotoz=lib.set_clphotoz_priors_mpp
setprior_clusteringphotoz.argtypes =[Double10,Double10]

get_sigma_8 = lib.get_sigma_8
get_h0 = lib.get_h0

log_like_wrapper = lib.log_like_wrapper

get_N_data = lib.get_N_data_masked
get_N_data.argtypes = []
get_N_data.restype = ctypes.c_int

get_N_tomo_shear = lib.get_N_tomo_shear
get_N_tomo_shear.argtypes = []
get_N_tomo_shear.restype = ctypes.c_int

get_N_tomo_clustering = lib.get_N_tomo_clustering
get_N_tomo_clustering.argtypes = []
get_N_tomo_clustering.restype = ctypes.c_int

get_N_ggl = lib.get_N_ggl
get_N_ggl.argtypes = []
get_N_ggl.restype = ctypes.c_int

initcmb = lib.init_cmb
initcmb.argtypes = [ctypes.c_char_p]

class IterableStruct(ctypes.Structure):
    def names(self):
        out = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                out.append(name)
            else:
                for i in range(length):
                    out.append(name + "_" + str(i))
        return out


    def iter_parameters(self):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                yield name, obj, ptype._length_
            else:
                yield name, obj, 0

    def iter_parameters_filter(self, used):
        for (name, obj, length) in self.iter_parameters():
            if name in used:
                yield name, obj, 0


    def convert_to_vector(self):
        p = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                p.append(obj)
            else:
                for i in range(length):
                    p.append(obj[i])
        return p

    def convert_to_vector_filter(self, used):
        p = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                if name in used:
                    p.append(obj)
            else:
                for i in range(length):
                    if name+'_'+str(i) in used:
                        p.append(obj[i])
        return p



    def read_from_cosmosis(self, block):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                for i in range(ptype._length_):
                    obj[i] = block[self.section_name, name+"_"+str(i)]
            else:
                setattr(self, name, block[self.section_name, name])



    def print_struct(self):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                for i in range(ptype._length_):
                    print ("%s[%d] = %f" % (name, i, obj[i]))
            else:
                print ("%s = %e" % (name, obj))


    def number_of_doubles(self):
        n=0
        for name, ptype in self._fields_:
            if hasattr(ptype, "_length_"):
                n += ptype._length_
            else:
                n += 1
        return n

    def set_from_vector(self, p):
        i=0
        j=0
        while i<len(p):
            name,ptype = self._fields_[j]
            j+=1
            if ptype == double:
                setattr(self, name, p[i])
                i+=1
            else:
                x = getattr(self, name)
                assert x._type_==double
                for k in range(x._length_):
                    x[k] = p[i]
                    i+=1


class InputCosmologyParams(IterableStruct):
    section_name = "cosmological_parameters"
    _fields_ = [
        ("omega_m", double),
        ("sigma_8", double),
        ("A_s", double),
        ("n_s", double),
        ("w0", double),
        ("wa", double),
        ("omega_b", double),
        ("omega_nuh2", double),
        ("h0", double),
        ("MGSigma", double),
        ("MGmu", double),
        ("theta_s", double),
    ]

    @classmethod
    def fiducial(cls):
        c = cls()
        c.omega_m = -1.0
        c.sigma_8 = -1.0
        c.A_s = 0.0
        c.n_s = 0.0
        c.w0 = 0.0
        c.wa = 0.0
        c.omega_b = 0.0
        c.omega_nuh2 = 0.0
        c.h0 = 0.0
        c.MGSigma = 0.0
        c.MGmu = 0.0
        c.theta_s = 0.0
        return c

    @classmethod
    def fiducial_sigma(cls):
        c = cls()
        c.omega_m = 0.02
        c.sigma_8 = 0.04
        c.A_s = 2.e-10
        c.n_s = 0.01
        c.w0 = .1
        c.wa = 0.02
        c.omega_b = 0.001
        c.omega_nuh2 = 0.0001
        c.h0 = 0.01
        c.MGSigma = 0.1
        c.MGmu = 0.1
        c.theta_s = 0.01
        return c




class InputNuisanceParams(IterableStruct):
    section_name = "nuisance_parameters"
    _fields_ = [
        ("bias", double*10),
        ("b_mag", double*10),
        ("lens_z_bias", double*10),
        ("source_z_bias", double*10),
        ("shear_m", double*10),
        ("p_ia", double*10)
    ]
    @classmethod
    def fiducial(cls):
        c = cls()
        c.bias[:] = [1.7, 1.7, 1.7, 2.0,2.0,2.0,2.0,2.0,2.0,2.0]
        c.b_mag[:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        c.lens_z_bias[:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        c.source_z_bias[:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        c.shear_m[:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        c.p_ia[:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        return c

    @classmethod
    def fiducial_sigma(cls):
        c = cls()
        c.bias[:] = np.repeat(0.2, 10)
        # c.bias2[:] = np.repeat(0.05, 10)
        c.b_mag[:] = np.repeat(0.05, 10)
        # c.b_ta[:] = np.repeat(0.05, 10)
        c.lens_z_bias[:] = np.repeat(0.002, 10)
        c.source_z_bias[:] = np.repeat(0.002, 10)
        c.shear_m[:] = np.repeat(0.002, 10)
        c.p_ia[:] = np.repeat(0.1, 10)
        return c

write_cosmolike_datavector.argtypes = [ctypes.c_char_p,InputCosmologyParams, InputNuisanceParams]
lib.log_like_wrapper.argtypes = [InputCosmologyParams, InputNuisanceParams]
lib.log_like_wrapper.restype = double
lib.get_sigma_8.argtypes=[InputCosmologyParams]
lib.get_sigma_8.restype = double

lib.get_h0.argtypes=[InputCosmologyParams]
lib.get_h0.restype = double

class LikelihoodFunctionWrapper(object):
    def __init__(self, varied_parameters, cosmo_min, cosmo_fiducial, cosmo_max, nuisance_min, nuisance_fiducial, nuisance_max):
        self.varied_parameters = varied_parameters
        self.cosmo_min = cosmo_min
        self.cosmo_fid = cosmo_fiducial
        self.cosmo_max = cosmo_max
        self.nuisance_min = nuisance_min
        self.nuisance_fid = nuisance_fiducial
        self.nuisance_max = nuisance_max



    def fill_varied(self, icp, inp, x):
        assert len(x) == len(self.varied_parameters), "Wrong number of parameters"
        i = 0
        for s in [icp, inp]:
            for name, obj, length in s.iter_parameters():
                if length==0:
                    if name in self.varied_parameters:
                        setattr(s, name, x[i])
                        i+=1
                else:
                    for j in range(length):
                        name_i = name + "_" + str(j)
                        if name_i in self.varied_parameters:
                            obj[j] = x[i]
                            i+=1

    def prior_cosmology(self, InputCosmologyParams):
        good = True
        for p in InputCosmologyParams.names():
            if p in self.varied_parameters:
                v = getattr(InputCosmologyParams,p)
                min_v = getattr(self.cosmo_min, p)
                max_v = getattr(self.cosmo_max,p)
                if v<min_v or v>max_v:
                #    print "Cosmo param {} out of bounds {}   [{},{}]".format(p,v,min_v,max_v)
                    good=False
        if good:
            return 0.0
        else:
            return -np.inf

    def prior_nuisance(self, InputNuisanceParams):
        params = ["bias", "b_mag", "p_ia"]
        good = True
        for p in params:
            for i in range(10):
                if '%s_%d'%(p,i) in self.varied_parameters:
                    min_val = getattr(self.nuisance_min,p)[i]
                    max_val = getattr(self.nuisance_max,p)[i]
                    value = getattr(InputNuisanceParams,p)[i]
                    if value<min_val or value>max_val:
                     #   print "Nuisance parameter {}[{}] out of bounds {}  [{},{}]".format(p,i,value,min_val, max_val)
                        good=False
        if good:
            return 0.0
        else:
            return -np.inf



    def __call__(self, x):
        icp = copy.deepcopy(self.cosmo_fid)
        inp = copy.deepcopy(self.nuisance_fid)
        self.fill_varied(icp, inp, x)
#        print
#        icp.print_struct()
#        inp.print_struct()
#        print
        flat_prior = self.prior_cosmology(icp) + self.prior_nuisance(inp)
        if not np.isfinite(flat_prior):
        #    print "outside flat prior range"
            return -np.inf,-1.
        like = lib.log_like_wrapper(icp, inp)
        print ("Likelihood = {}".format(like))
        # print inp.p_ia[0], inp.p_ia[1], '<-- p_ia'
        # print inp.bias[0], inp.bias[1], inp.bias[2], inp.bias[3], inp.bias[4], '<-- bias'
        # print inp.b_mag[0], inp.b_mag[1], inp.b_mag[2], inp.b_mag[3], inp.b_mag[4], '<-- bmag'
        if like < -1.0e+14:
            return -np.inf,-1.

        return like,get_sigma_8(icp)

# def sample_cosmology_only(MG = False):
#     if MG:
#         varied_parameters = InputCosmologyParams().names()
#     else:
#         varied_parameters = ['omega_m']
#         varied_parameters.append('sigma_8')
#         varied_parameters.append('n_s')
#         varied_parameters.append('w0')
#         varied_parameters.append('wa')
#         varied_parameters.append('omega_b')
#         varied_parameters.append('h0')
#         varied_parameters.append('theta_s')

#     return varied_parameters

# def sample_cosmology_6x2_allsys(tomo_N_shear,tomo_N_lens,MG = False):
#     varied_parameters = sample_cosmology_only(MG)
#     varied_parameters += ['bias_%d'%i for i in xrange(tomo_N_lens)]
#     varied_parameters += ['b_mag_%d'%i for i in xrange(tomo_N_lens)]
#     varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
#     # varied_parameters.append('source_z_s')
#     varied_parameters += ['lens_z_bias_%d'%i for i in xrange(tomo_N_lens)]
#     # varied_parameters.append('lens_z_s')
#     varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
#     # 2 parameters for IA NLA model
#     varied_parameters += ['p_ia_%d'%i for i in xrange(2)]


#     # varied_parameters.append('A_ia')
#     # varied_parameters.append('beta_ia')
#     # varied_parameters.append('eta_ia')
#     # varied_parameters.append('eta_ia_highz')
#     # varied_parameters += ['bary_%d'%i for i in xrange(3)]
#     return varied_parameters



#def likelihood_task(p):
#    return my_likelihood(p)

def sample_main(varied_parameters, iterations, nwalker, nthreads,
        filename, cosmo_min, cosmo_fid, cosmo_max, nuisance_min, nuisance_fid, nuisance_max,
        pool=None):
    cosmo_fid.print_struct()
    nuisance_fid.print_struct()
    print (varied_parameters)

    likelihood = LikelihoodFunctionWrapper(varied_parameters,
        cosmo_min, cosmo_fid, cosmo_max, nuisance_min,
        nuisance_fid,nuisance_max)
#    global my_likelihood
#    my_likelihood = likelihood
    starting_point = []
    starting_point += cosmo_fid.convert_to_vector_filter(varied_parameters)
    starting_point += nuisance_fid.convert_to_vector_filter(varied_parameters)

    print ("Starting point = ", starting_point)

    std = InputCosmologyParams.fiducial_sigma().convert_to_vector_filter(varied_parameters)
    std += InputNuisanceParams().fiducial_sigma().convert_to_vector_filter(varied_parameters)

    p0 = emcee.utils.sample_ball(starting_point, std, size=nwalker)

    ndim = len(starting_point)
    print ("ndim = ", ndim)
    print ("start = ", starting_point)
    print ("std = ", std)
    #dtype = [("h", float), ("sigma_8", float)]
    sampler = emcee.EnsembleSampler(nwalker, ndim, likelihood, #blobs_dtype=dtype,
        threads=nthreads, pool=pool)
    ### try to write to multiple files
    #try:
    #    from mpi4py import MPI
    #    pid = MPI.COMM_WORLD.rank
    #    psize = MPI.COMM_WORLD.size
    #    filename = filename+"_%d-%d"%(pid, psize)
    #except ImportError:
    #    print("mpi4py not found\n Write chains into single file\n")
    #f = open(filename, 'w')
    f = open(filename, 'w')

    #write header here
    f.write('# ' + '    '.join(varied_parameters)+"  sigma_8  log_like\n")


    for (p, loglike, state, blobs) in sampler.sample(p0,iterations=iterations):
#        for row,h,s8,logl in zip(p,blobs["h"],blobs["sigma_8"],loglike):
        for row,s8,logl in zip(p,blobs,loglike):
            p_text = '  '.join(str(r) for r in row)
            f.write('%s  %e %e\n' % (p_text,s8,logl))
        f.flush()
    f.close()

    pool.close()
