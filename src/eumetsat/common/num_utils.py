'''
Created on Feb 14, 2011

@author: guillaume.aubert@eumetsat.int
'''
import struct
import string

def hex2bin(str):
   bin = ['0000','0001','0010','0011',
         '0100','0101','0110','0111',
         '1000','1001','1010','1011',
         '1100','1101','1110','1111']
   aa = ''
   for i in range(len(str)):
      aa += bin[string.atoi(str[i],base=16)]
   return aa

def float_to_bin(a_float):
    """ convert int into a 32 bits binary string representation """
    
    s = struct.pack('>f',a_float)
    hexa = ''.join('%2.2x' % ord(c) for c in s)
    
    return hex2bin(hexa)
    
def float_to_hex(a_float):
    """convert a float in hexa """
    s = struct.pack('>f',a_float)
    hexa = ''.join('%2.2x' % ord(c) for c in s)
    return hex2bin(hexa)



if __name__ == '__main__':
                                     
    assert float_to_bin(0.15625)  == "00111110001000000000000000000000"
    
    assert float_to_bin(-0.15625) == "10111110001000000000000000000000"
    
    assert float_to_bin(173.3125) == "01000011001011010101000000000000"