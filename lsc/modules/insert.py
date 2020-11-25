#!/usr/local/bin/python
#


list_directory = "/p/eaddata1/lpku/TSC/v10.8f/lsc/list" 
in_directory = "/p/eaddata1/lpku/TSC/v10.8f/lsc/modules/WORK" 
out_directory = "/p/eaddata1/lpku/TSC/v10.8f/lsc/modules" 
list_file_name = list_directory + "/" + "module_list"

f = open(list_file_name, "r")

sub_list = []

while 1:
    line = f.readline()
    if (not line) : break

    sub_list.append(line[0:-1])

f.close()


for file in sub_list:

    print 'processing file: ', file

    in_file = in_directory + "/" + file
    out_file = out_directory + "/" + file
    fi = open(in_file, "r")
    fo = open(out_file, "w")

    module_name = file[file.find('_')+1:file.find('.')]

    lower = 'abcdefghijklmnopqrstuvwxyz0123456789'
    upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'

    new_name = ""
#   for i in range(len(module_name)) :
#      new_name = new_name + upper[lower.find(module_name[i])]
    new_name = module_name
	
    lines = []
    lines = fi.readlines()
    no_of_lines = len(lines)

    s = '      MODULE ' + new_name +"\n"
    fo.write(s)
    s = '      USE PARAMS' +"\n"
    fo.write(s)
    s = '      USE EMPARAMS' +"\n"
    fo.write(s)
    s = '      IMPLICIT NONE' +"\n"
    fo.write(s)
 
    find = 0
    for s in lines : 
        if s.find("COMMON") is -1 :
            if find : 
                if len(s) < 5 : 
                    fo.write(s)
                else :
                    if s[5] is not " " :
                        ss = "!" + s[1:]
                        s = ss
                        fo.write(s)
                    else :
                        find = 0
                        fo.write(s)
            else :
                find = 0
                fo.write(s)
        else :
            ss = "!" + s[1:]
            find = 1
            s = ss
            fo.write(s)

    s = '      END MODULE ' + new_name +"\n"
    fo.write(s)


    fi.close()
    fo.close()

# finish
        
            
    
