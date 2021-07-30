from collections import defaultdict
from dwave.system import LazyFixedEmbeddingComposite, EmbeddingComposite, DWaveSampler
import dwave.inspector
import pdb
import scipy.io
import numpy


# J can be different per channel use Matlab script version 13

mat = scipy.io.loadmat('inVars2.mat')

ag=len(mat['Gt2'])
bg=len(mat['Gt2'][0])
eg=len(mat['Gt2'][0][0])
cg=len(mat['Gt2'][0][0][0])
try:
    dg=len(mat['Gt'][0][0][0][0])
except:
    dg=1
    probdim = 3

af=len(mat['Ft'])           #number of variables
bf=len(mat['Ft'][0])        #frame size
cf=len(mat['Ft'][0][0])     #number of frames
try:
    df=len(mat['Ft'][0][0][0])  #snr points
except:
    df=1
    probdim = 3

flagIter = 0 #0 on first iteration

#bf=10 ##################

for ss in range(df):
    for ff in range(cf):

        print("Frame ", ff+1, "\\", cf, " SNRpoint ", ss+1, "\\", df, "\n")
        #pdb.set_trace()

        #G = [ [0]*bg for i in range(ag)]
        #for aa in range(ag):
        #    for bb in range(bg):
        #        if probdim == 3: #single SNR point
        #            G[aa][bb] = mat['Gt'][aa][bb][ff]
        #        else:
        #            G[aa][bb] = mat['Gt'][aa][bb][ff][ss]

        G2 = [ [ [0]*eg for i in range(ag)] for j in range(bg)]
        for ee in range(eg):
            for aa in range(ag):
                for bb in range(bg):
                    if probdim == 3: #single SNR point
                        G2[aa][bb][ee] = mat['Gt2'][aa][bb][ee][ff]
                    else:
                        G2[aa][bb][ee] = mat['Gt2'][aa][bb][ee][ff][ss]

        F = [ [0]*bf for i in range(af)]
        for aa in range(af):
            for bb in range(bf):
                if probdim ==3:
                    F[aa][bb] = mat['Ft'][aa][bb][ff]
                else:
                    F[aa][bb] = mat['Ft'][aa][bb][ff][ss]

   

        h = defaultdict(float)
        J = defaultdict(float)
        d = af
        blk = bf
        #pdb.set_trace()
        #for bb in range(blk):
        #    for r in range(d):
        #        for c in range(d):
        #            J[(bb*d+r,bb*d+c)] = G[r][c]

        #for bb in range(blk):
        #    for r in range(d-1):
        #        for c in range(r+1,d):
        #            J[(bb*d+r,bb*d+c)] = G[r][c]

        for bb in range(blk):
            for r in range(d-1):
                for c in range(r+1,d):
                    J[(bb*d+r,bb*d+c)] = G2[r][c][bb]

        vec = bf
        for vv in range(vec):
            for r in range(d):
                h[vv*d+r] = F[r][vv]

        #pdb.set_trace()   

        if flagIter == 0:
            # Define the sampler that will be used to run the problem
            dsampler =  DWaveSampler(solver={'qpu':True, 'topology__type': 'pegasus'})
            #sampler = EmbeddingComposite(dsampler)
            sampler = LazyFixedEmbeddingComposite(dsampler)
            #sampler = EmbeddingComposite(DWaveSampler(solver={'qpu':True}))

        #pdb.set_trace()

        # Run the problem on the sampler and print the results
        sampleset = sampler.sample_ising(h, J,
                                        num_reads = 50,
                                        annealing_time = 20,
                                        reduce_intersample_correlation = False,
                                        label='8x8_16qam')
        #print(type(sampleset))
        #print(sampleset)

        rslt = sampleset.record
        if flagIter == 0:
            rsltAcc = rslt
        else:
            rsltAcc = numpy.hstack((rsltAcc,rslt))
        flagIter = 1
        numpy.save('outVarsT.npy', rsltAcc)

        dwave.inspector.show(sampleset)
        pdb.set_trace()
        
scipy.io.savemat('outVarsT.mat', mdict={'rsltAcc': rsltAcc})
dwave.inspector.show(sampleset)
pdb.set_trace()