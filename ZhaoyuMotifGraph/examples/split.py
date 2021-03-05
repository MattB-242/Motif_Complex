import re
import linecache

def split(input_file_path,output_file_directory,keyword='END',outfile_name='Crystal'):
    fp = open(input_file_path, 'r')

    number = []
    lineNumber = 1

    for eachLine in fp:
        m = re.search(keyword, eachLine)  ##查询关键字
        if m is not None:
            number.append(lineNumber)  # 将关键字的行号记录在number中
        lineNumber = lineNumber + 1
    size = int(len(number))

    start = 1
    end = number[0]
    destLines = linecache.getlines(input_file_path)[start:end - 1]
    fp_w = open(output_file_directory + outfile_name + '_' +0 + '.cif', 'w')
    for key in destLines:
        fp_w.write(key)
    fp_w.close()

    for i in range(0, size - 1):
        start = number[i]
        end = number[i + 1]
        destLines = linecache.getlines(input_file_path)[start + 1:end - 1]
        fp_w = open(output_file_directory + outfile_name + '_' + str(i + 1) + '.cif', 'w')
        for key in destLines:
            fp_w.write(key)
        fp_w.close()


