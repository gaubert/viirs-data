'''
Created on Mar 8, 2011

@author: guillaume.aubert@eumetsat.int
'''


# Initialize the module
import eugene
import datetime
import h5py

A_DATE = datetime.datetime.now()



def extract_record_content(record_name, hdf_group):
    """ generically extract a record content """
    
    # create record
    record = prod.get(record_name)
    group  = hdf_file.create_group(hdf_group)
    
    for (index,key) in enumerate(record._children_):
        if key not in ("RECORD_HEADER"):
            # forced to do a eval because the API doesn't provide anything to iter over the record attributes: ggrrrrr
            val = eval('record.%s' %(key))
            
            print("key = %s, val = %s, type(val) = %s\n" %( key, val, type(val)))
            
            if type(val) == type(A_DATE):
                group.create_dataset(key, data=str(val))
            else:
                group.create_dataset(key, data=val)


if __name__ == '__main__':
    
    dir = "/homespace/gaubert/eugene-4.9.2/inst/bin"
    filename = "GOME_earthshine_xxx_1B_M02_20101124143858Z_20101124144158Z_N_O_20101124161022Z"
    
    # Load a product file
    prod = eugene.Product("%s/%s" % (dir, filename))
    
    # create hdf5 file
    hdf_file = h5py.File('/tmp/myfile.hdf5','w')
    
    for rec in prod._records_:
        print("record = %s" % (rec))
        if rec not in ('mdr-1b-earthshine', 'mdr-1b-moon', 'mdr-1b-sun', 'mdr-1b-calibration', 'VEADR-1', 'viadr-smr', 'dummy', 'IPRS'):
            extract_record_content(rec,rec)
        else:
            print("IGNORE %s\n" %(rec))
    
    #extract mphr
    #extract_record_content("MPHR","MPHR")
    
    #extract  sphr
    #extract_record_content("SPHR","SPHR")
    
    #extract_record_content("geadr-elevation", "geadr-elevation")
   