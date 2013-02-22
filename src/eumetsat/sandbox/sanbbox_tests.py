'''
Created on Mar 11, 2011

@author: gaubert
'''

import eumetsat.common.num_utils as num_utils

def print_decomposed_values(val):
    """ print properly the values """
    #print num_utils.float_to_bin(-21.1949)
    decomp =  num_utils.decomposed_rep_float_to_bin(val)
    print("%s %s %s : %f" %(decomp['sign'], decomp['exponent'], decomp['mantissa'], val))

if __name__ == '__main__':
    
    print("s  exp      mantissa")
    val1 = -21.1949
    print_decomposed_values(val1)
    
    val2 = -21.195
    print_decomposed_values(val2)
    
    print()
    
    print_decomposed_values(1.0)
    
    val3 = -21.2576
    print_decomposed_values(val3)
    
    val4 = -21.258
    print_decomposed_values(val4)
    
    
    print("Orig following values\n")
    print_decomposed_values(-21.1949)
    print_decomposed_values(-21.2076)
    print_decomposed_values(-21.2202)
    print_decomposed_values(-21.2328)
    print_decomposed_values(-21.2452)
    print_decomposed_values(-21.2576)
    
    print("rounded values")
    print_decomposed_values(-21.195)
    print_decomposed_values(-21.208)
    print_decomposed_values(-21.22)
    print_decomposed_values(-21.233)
    print_decomposed_values(-21.245)
    print_decomposed_values(-21.258)
    
    
    
    
    