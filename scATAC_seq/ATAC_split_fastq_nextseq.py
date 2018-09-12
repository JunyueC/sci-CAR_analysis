# python ATAC_split_fastq.py read1 read2 output_folder P7_index P5_index
'''
This script accept a read1, read2 file, a barcode file P5 and P7 file, output
file folder, then it read in the read1 and read2 file, split the read1 and read2 based on
the combination of the P5 and P7 barcode, and attach the adaptor barcode to the read1 and read2 
sequence
'''

'''
First, I am going to create the files with the P5-P7 index combination,

For each read1, I will go through line by line, for line with the index line,
extract the index, and file which this index belong to, and then write the read1 and read2 both into the 
output file.
'''
import subprocess
import sys
from Levenshtein import distance
import gzip

def extract_barcode(index):
    '''
    This function accept a index sequence, and return the extracted P7+P5 barcode
    in the index, and the barcode in the Tn5 adaptor
    '''
    barcodes = index.split('+')
    P7 = barcodes[0][8:]
    P5 = barcodes[1][8:]
    N7 = barcodes[0][0:8]
    N5 = barcodes[1][0:8]
    
    P_index = P7 + '.' + P5
    N_index = N7 + '.' + N5
    
    return (P_index, N_index)

'''
def P7_P5_list(P7_file, P5_file):
    
    #this function accept a text including P7 index and a text file including a 
    #P5 index and then generate a list for all the P7.P5 barcode combinations
    
    P7 = open(P7_file)
    P5 = open(P5_file)
    index1 = P7.readlines()
    index2 = P5.readlines()
    index1 = map(lambda x: x.strip(), index1)
    index2 = map(lambda x: x.strip(), index2)
    all_indexes = []
    for i1 in index1:
        for i2 in index2:
            all_indexes.append(i1 + '.' + i2)
    
    P7.close()
    P5.close()
    return all_indexes
'''

def P7_P5_list(P7_file, P5_file):
    '''
    this function accept a text including P7 index and a text file including a 
    P5 index and then generate a list for all the P7.P5 barcode combinations
    '''
    P7 = open(P7_file)
    P5 = open(P5_file)
    index1 = P7.readlines()
    index2 = P5.readlines()
    index1 = map(lambda x: x.strip(), index1)
    index2 = map(lambda x: x.strip(), index2)
    all_indexes = []
    for i in range(len(index1)):
            all_indexes.append(index1[i] + '.' + index2[i])
    
    P7.close()
    P5.close()
    return all_indexes


def find_P7_P5_barcode(index, P7_P5_barcodes):
    '''
    this function accept a P7_P5 barcode in the read, and a list of P7_P5 barcodes provided,
    and then find the barcode in the P7_P5 barcode list that match the index with distance <2,
    then return the found barcode; if no barcode is found, then return -1
    '''
    result = -1

    for barcode in P7_P5_barcodes:
        diff = distance(index, barcode)
        if (diff <= 1):
            result = barcode
            break
    
    return result
    
def split_fastq(read1, read2, output_folder, P7_P5_barcodes): 
    '''
    This function accept a read1 file, a read2 file, a list of P7-P5 barcodes, create the output files,
    and then go through each line of the read1 and read2 file,
    and then check the barcode in the read1, compare it with the P7_P5 barcodes,
    if the distance < 2, then output the reads into the output files
    '''
    
    # generate the output files
    output_read1 = {}
    output_read2 = {}
    sample_ID = open(output_folder + "/sample_ID.txt", 'w')
    barcode_count = {}
    for barcode in P7_P5_barcodes:
        output_R1 = output_folder + '/' + barcode + '.R1.fastq.gz'
        output_R2 = output_folder + '/' + barcode + '.R2.fastq.gz'
        output_read1[barcode] = gzip.open(output_R1, 'w')
        output_read2[barcode] = gzip.open(output_R2, 'w')
        sample_ID.write(barcode + '\n')
        barcode_count[barcode] = 0
        
    sample_ID.close()
    # open the read1 and read2 file, and then for each read, check the barcode and output the matched read into the
    # output file
    
    f1 = gzip.open(read1, 'rb')
    f2 = gzip.open(read2, 'rb')
    line_n = 0
    total_read=0
    
    while(True):
        line1 = f1.readline()
        line2 = f2.readline()
        line_n += 1
        if (not line1):
            break
        if (line_n % 4 == 1):
            index = line1.strip().split(':')[-1]
            #print "index: ", index
            index = extract_barcode(index)
            P_index = index[0]
            N_index = index[1]
            P_index = find_P7_P5_barcode(P_index, P7_P5_barcodes)
            total_read+=1
            if (P_index != -1):
                output_R1 = output_read1[P_index]
                output_R2 = output_read2[P_index]
                barcode_count[P_index] += 1
                line1 = '@' + N_index + ',' + line1[1:]
                line2 = '@' + N_index + ',' + line2[1:]
                output_R1.write(line1)
                output_R2.write(line2)
                
                # output the other lines into the files
                line1 = f1.readline()
                line2 = f2.readline()
                output_R1.write(line1)
                output_R2.write(line2)
                line_n += 1
                
                line1 = f1.readline()
                line2 = f2.readline()
                output_R1.write(line1)
                output_R2.write(line2)
                line_n += 1
                
                line1 = f1.readline()
                line2 = f2.readline()
                output_R1.write(line1)
                output_R2.write(line2)
                line_n += 1
                
            else:
                line1 = f1.readline()
                line2 = f2.readline()
                line_n += 1
                
                line1 = f1.readline()
                line2 = f2.readline()
                line_n += 1
                
                line1 = f1.readline()
                line2 = f2.readline()
                line_n += 1
    # close the files
    f1.close()
    f2.close()
    print "Total read count: ", total_read
    asigned_read = 0
    for barcode in P7_P5_barcodes:
        output_read1[barcode].close()
        output_read2[barcode].close()
        print "Read_count: ", barcode, barcode_count[barcode]
        asigned_read += barcode_count[barcode]
    print "Total asigned read: ", asigned_read

def split_fastq_main(read1, read2, output_folder, P7_file, P5_file):
    '''
    this is the main function that accept the read1 file, the read2 file, the output folder,
    the files including P7 barcode, P5 barcode and then split the 
    the read1 and read2 files based on the P5 and P7 combinations
    '''
    print "Start splitting the fastq files..."
    P7_P5_barcodes = P7_P5_list(P7_file, P5_file)
    print "Number of output files: ", len(P7_P5_barcodes)
    
    split_fastq(read1, read2, output_folder, P7_P5_barcodes)
    print "fastq files splitted~"

if __name__ == "__main__":
    read1 = sys.argv[1]
    read2 = sys.argv[2]
    output_folder = sys.argv[3]
    P7_file = sys.argv[4]
    P5_file = sys.argv[5]
    split_fastq_main(read1, read2, output_folder, P7_file, P5_file)