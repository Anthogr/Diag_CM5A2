def readncfile(inpath, indic, mydata):
    from netCDF4 import Dataset
    
    my_data=mydata
    f = Dataset(inpath)
    
    #ANTHOG added: Case insensitive variable name search
    #-----------------------------------------------#
    # For exemple if you entrered in a specific dictionnary a variable name like 
    # x = 'nav_lon' and it is actually 'NAV_LON' in the nc file the few line below 
    # will change x = 'NAV_LON' so the variable will be found and loaded
    
    dictKeys     = f.variables.keys()
    dictKeysList = list(dictKeys)
    
    mapVar            = (map(lambda x: x.casefold(), dictKeysList))
    dictKeysListLower = list(mapVar)
    #-----------------------------------------------#
    
    for x in indic:
        
        varname = x
        var2read = indic[x]
        
        if var2read in f.variables:
            
            my_data[varname] = f.variables[var2read][:]
            
        elif var2read.casefold() in dictKeysListLower: # ANTHOG Case insensitive variable name search condition
            
            ind      = dictKeysListLower.index(var2read.casefold())
            var2read = dictKeysList[ind]
            
            my_data[varname] = f.variables[var2read][:]
        
    f.close()
    return(my_data)