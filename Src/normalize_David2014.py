#!/usr/bin/env python 

import sys, getopt
import numpy as np
from math import log10, floor, pow

# Parsing parameters modified from http://www.tutorialspoint.com/python/python_command_line_arguments.htm
inputFile = ''
outputFile = ''
try:
	opts, args = getopt.getopt(sys.argv[1:],"i:o:")
except getopt.GetoptError:
	print 'normalize_David2014.py -i <inputFile.tsv> -o <outputFile.tsv>'
	# example: normalize_David2014.py -i 'otuTablePrenormDavid14.tsv' -o 'otuTablePostnormDavid14.tsv'
	sys.exit(2)
for opt, arg in opts:
	if opt == "-i":
		inputFile = arg
	elif opt == "-o": 
		outputFile = arg
## End form parsing

otuNp=np.loadtxt(inputFile, dtype='f', comments='#', delimiter='\t', converters=None, skiprows=0, usecols=None, unpack=False, ndmin=2)
time_c = len(otuNp) # num of samples
otu_c = len(otuNp[1]) # num of OTUs
seq_M = otuNp
############# ME until here



################# BEGIN code from http://nbviewer.ipython.org/github/ladavid/mit_timeseries/blob/master/NormalizeDemo.ipynb ##################
# Function definitions
def JSD(v1,v2): # Defined in https://github.com/ladavid/mit_timeseries/blob/master/ExampleHelper.py
    ''' compute jensen-shannon divergence '''
    half=(v1+v2)/2
    kl_1_v = v1*np.log2(v1/half)
    kl_2_v = v2*np.log2(v2/half)
    kl_1_v[np.isnan(kl_1_v)] = 0.
    kl_2_v[np.isnan(kl_2_v)] = 0.
    return 0.5*sum(kl_1_v)+0.5*sum(kl_2_v)
#enddef

def WeightedMedian(data_M,weight_M,i): # Defined in https://github.com/ladavid/mit_timeseries/blob/master/ExampleHelper.py
    ''' compute weighted median on data matrix.  modeled on
    R weighted.median function. '''
    # find samples
    sample_v = np.nonzero(~np.isnan(data_M[:,0]))[0]
    sample_i = np.nonzero(sample_v==i)[0][0]
    sample_M = data_M[sample_v,:]
    weight_v = weight_M[i,sample_v]
	
    # sort each matrix column
    sort_M = np.argsort(sample_M,0)
	
    # find halfway points for weights in each column
    otu_n = sort_M.shape[1]
    exp_v = np.zeros((1,otu_n))[0]
    for col_i in range(otu_n):
        mid_i = np.nonzero(np.cumsum(weight_v[sort_M[:,col_i]]) > 0.5)[0][0]
        median_n = sample_M[sort_M[mid_i,col_i],col_i]
        exp_v[col_i] = median_n
    #endfor
	    
    # return vector of values
    return exp_v
#enddef


# Compute sample pairwise similarity matrices using Jensen-Shannon Distance
import scipy.spatial.distance as ssd
seq_sum_v = np.sum(seq_M,1)  # sum all reads per sample (i.e. sum all columns per each row)
seq_frac_M = seq_M/np.tile(seq_sum_v,(otu_c,1)).T  # compute relative frequency, as TSS, i.e. dividing each absolute frequency (in seq_M) by the total number of reads per sample (in seq_sum_v). With np.tile duplicated the denominator (the total number of reads per sample) all the times is needed, it means, one per OTU.
# ME: Applying pseudocount: "a value smaller than the minimum abundance value before transformation." [https://www.researchgate.net/publication/261220701_A_fair_comparison, Costea,2014, NatureMethods, DOI: 10.1038/nmeth.2897, suppl.info: http://www.nature.com/nmeth/journal/v11/n4/extref/nmeth.2897-S1.pdf]"
minAbund = seq_frac_M[np.nonzero(seq_frac_M)].min() # >1.6233187e-06. Get the order of magnitude based on 10 bases, and take the next smaller number than this one. It means, in the example, for an order of magnitude of the minimum abundance of 1e-06, the pseudo-count will be 1e-07.
print 'min. Abundance=',minAbund
orderOfMagnitude=int(floor(log10(minAbund)))
# ME: To compute pseudocount automatically
pseudoCount = pow(10,orderOfMagnitude-1)
print 'pseudoCount=', pseudoCount
# ME: Sum pseudocount to ALL values in matrix
for row_i in range(time_c):
	for col_i in range(otu_c):
		seq_frac_M[row_i][col_i] = seq_frac_M[row_i][col_i] + pseudoCount
	#endfor col
#endfor row

#jsd_M = ssd.squareform(ssd.pdist(seq_frac_M,ExampleHelper.JSD)) # I copy the function here, outside the module 'ExampleHelper'
jsd_M = ssd.squareform(ssd.pdist(seq_frac_M,JSD)) # Compute the JSD distance and convert to an array with the bottom diagonal matrix, since the distance matrix is square and symmetric.
dist_M = np.sqrt(jsd_M) # Compute sqrt of the JSD distance


# Perform normalization
# Setup weighting matrix using JSD. Use simple scheme based on squares.
weight_M = np.zeros((time_c,time_c)) + np.nan
time_v = range(time_c)
for time_i in time_v:
    weight_v = np.power(1 - jsd_M[time_i,:],2)
    weight_v = weight_v/np.sum(weight_v)
    weight_M[time_i,:] = weight_v
    weight_M[time_i,time_i] = 0.
#endfor

# Simplify by looking at most abundant OTUs (defined as the fewest number tha account for 90% of reads)
#In [12]:
median_v = np.median(seq_M,0)
frac_v = median_v/sum(median_v)
sort_v = np.flipud(np.sort(frac_v))
cumsum_v = np.cumsum(sort_v)
thresh_i = np.nonzero(cumsum_v > 0.9)[0][0]
otu_subset_v = np.nonzero(frac_v > sort_v[thresh_i])[0]

#Start the normalization using fractional abundances converted to log-space
#In [13]:
obs_M = np.log10(seq_frac_M[:,otu_subset_v])
start_M = obs_M.copy()

#Normalize each time point, in random order
#In [14]:
np.random.shuffle(time_v)
for t in time_v:
	# observed abundances
	obs_v = obs_M[t,:]
	        
	# expected abundances
	#exp_v = ExampleHelper.WeightedMedian(obs_M,weight_M,t)    
	exp_v = WeightedMedian(obs_M,weight_M,t)    
	
	# solve y=mx in log-space, using median, so as to avoid effects of outliers
	m_t = np.median(obs_v-exp_v)
	   
	# update time-series
	obs_M[t,:] = obs_M[t,:] - m_t    
#endfor


#Now, rescale simulated sequence data using inferred scaling factors.
#In [15]:
scale_v = start_M[:,0]-obs_M[:,0] # scaling factors for each time point
scale_M = np.tile(scale_v,(seq_frac_M.shape[1],1)).T
scale_frac_M = np.log10(seq_frac_M) - scale_M
inf_abun_M = np.power(10,scale_frac_M) # inferred OTU abundances in linear space
################# END code from http://nbviewer.ipython.org/github/ladavid/mit_timeseries/blob/master/NormalizeDemo.ipynb ##################

np.savetxt(outputFile, inf_abun_M, fmt='%.18e', delimiter='\t', newline='\n', header='', footer='', comments='# ')



