#!/usr/local/bin/python

list_directory = "/p/eaddata1/lpku/TSC/v10.8f/lsc/list" 
in_directory = "/p/eaddata1/lpku/TSC/v10.8f/lsc/modules" 
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

    subroutine_name = file.replace(".f",".f90")
    in_file = in_directory + "/" + file
    out_file = out_directory + "/" + subroutine_name
    fi = open(in_file, "r")
    fo = open(out_file, "w")


    lines = []
    lines = fi.readlines()
    no_of_lines = len(lines)

    lower = 'abcdefghijklmnopqrstuvwxyz0123456789'
    upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
    blank = " "
    amper = "&"
    sa = lines[0]
    sb = lines[1]

    line_count = 0

    while line_count < no_of_lines-1 :

        line_count = line_count + 1
        found = 0
        comment =""
        skip = 0

        line = lines[line_count]
        ls = line_count
        while ls < no_of_lines and (line[0] is not blank or len(line) < 5) :
            skip = skip + 1
            ls = ls + 1
            sb = ""
            if ls < no_of_lines : sb = lines[ls]
            line = sb

        sb = line
        if len(line) > 5 :
            if line[5] is not blank :
                found = 1
                sb = "     &" + line[6:]

        if found is 1 :
            line = sa
            send = line.find("!")
            if send is -1 :
                sa = line.replace("\n"," ")
                ex = 72-len(sa)
                ex = max(0,ex)
                line = sa[:len(sa)] + blank*ex + " &  \n"
                sa = line
            else :
                sa = line[:send]
                comment = blank*(send-1) + line[send:]
                ex = 72-len(sa)
                ex = max(0,ex)
                line = sa + blank*ex + " &  \n"
                sa = line

        line = sa
        if line[0] is "*" : sa = "!" + line[1:] 
        if line[0] is "c" : sa = "!" + line[1:] 
        if line[0] is "C" : sa = "!" + line[1:] 
        fo.write(sa)
        if len(comment) is not 0 : fo.write(comment)
        i = 0
        while i < skip :
            i = i + 1
            sa = lines[line_count]
            line = sa
            if line[0] is "*" : sa = "!" + line[1:]
            if line[0] is "c" : sa = "!" + line[1:]
            if line[0] is "C" : sa = "!" + line[1:]


            fo.write(sa)
            line_count = line_count + 1
        sa = sb

    fo.write(sa)
    fo.close()

# finish
        
            
    
