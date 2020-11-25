#!/usr/local/bin/python
#

in_directory = "/p/eaddata1/lpku/TSC/v10.8f/lsc/include" 
out_directory = "/p/eaddata1/lpku/TSC/v10.8f/lsc/modules" 
list_file_name = out_directory + "/" + "junk"

f = open(list_file_name, "r")

sub_list = []

while 1:
    line = f.readline()
    if (not line) : break

    sub_list.append(line[0:-1])

f.close()


for file in sub_list:

    print 'processing file: ', file

    subroutine_name = "module_" + file.replace(".inc",".f")
    in_file = in_directory + "/" + file
    out_file = out_directory + "/" + subroutine_name
    fi = open(in_file, "r")
    fo = open(out_file, "w")


    lower = 'abcdefghijklmnopqrstuvwxyz0123456789'
    upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'

    lines = []
    lines = fi.readlines()
    no_of_lines = len(lines)
    line_count = 0
    while line_count < no_of_lines :
        line = lines[line_count]
        fo.write(line)
        line_count = line_count + 1

    fi.close()
    fo.close()

# finish
        
            
    
