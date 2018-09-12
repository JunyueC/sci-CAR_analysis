"""
Created on Mon Aug 10 22:07:08 2015

@author: Junyue
"""

'''
python intersection_dnase.py input_folder sampleID Dnase_hot_spot output_folder
this script take the input folding containing bed file and sampleID and danse hot spot
bed file, and generate the start site for all the bed file and calculate the
percentage of intersection for all the bed file and output the result to the output_folder
'''

import sys
import subprocess
def intersection_dnase(input_folder, sampleID, Dnase_hot_spot, output_folder):
    print '''
    *******************calculate the intersection between bed file and Dnase hot spot **********************
    input folder: %s
    sampleID: %s
    Dnase hot spot: %s
    output folder: %s
    ********************************************************************************************************
    ''' %(input_folder, sampleID, Dnase_hot_spot, output_folder)
    
    #make output directory
    make_output_folder = 'mkdir -p ' + output_folder
    err = subprocess.call(make_output_folder, shell=True)
    
    #make output directory for the start site of the bed file
    start_site_folder = input_folder + '/' + 'start_site'
    make_start_site_folder = 'mkdir -p ' + start_site_folder
    err = subprocess.call(make_start_site_folder, shell=True)
    
    #loop through all the samples, generage the start site bed file and calculate
    # the intersection between the bed file and the Dnase hot spot
    f1 = open(sampleID, 'r')
    for i in f1:
        sample = i.strip()
        input_file = input_folder + '/' + sample + '.bed'
        start_site_file = start_site_folder + '/' + sample + '_start_site.bed'
        output_file = output_folder + '/' + sample + '_interesection_hot_spot.txt'
        make_start_site_file = '''awk '{if($6=="+"){print $1"\t"$2-1"\t"$2"\t"$4} else{print $1"\t"$3"\t"$3+1"\t"$4}}' ''' \
        + input_file + ' >' + start_site_file
        err = subprocess.call(make_start_site_file, shell=True)
        
        if err != 0:
            raise IOError("!!!!!!!!!!!!!!program error, cannot generate start site bed file!!!!!!!!!!!!!")
        else:
            print "~~~~~~~~~~~~~~sample: %s. start site bed file generated~~~~~~~~~~~~~" %(sample)
        
        
        make_intersection = 'echo ' + sample + ' >' + output_file + '; bedtools intersect -a ' + start_site_file + ' -b ' + Dnase_hot_spot \
        + ' -wa  -nonamecheck |sort -u - | wc -l - >>' + output_file + '; sort -u ' + start_site_file \
        + ' |wc -l - >>' + output_file + '; echo ----------------------- >>' + output_file
        
        print "Here is how to make the intersection: %s" % (make_intersection)
        err = subprocess.call(make_intersection, shell=True)
        
        if err != 0:
            raise IOError("!!!!!!!!!!!!program error. Intersection calculation failed!!!!!!!!!!!!!!!!!")
        else:
            print "~~~~~~~~~~~~~~~~~~~~sample: %s. interesction calculation done~~~~~~~~~~~~~~~~~~" %(sample)
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~all intersection calculation done~~~~~~~~~~~~~~~~~~~"

if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    Dnase_hot_spot = sys.argv[3]
    output_folder = sys.argv[4]
    intersection_dnase(input_folder, sampleID, Dnase_hot_spot, output_folder)
    