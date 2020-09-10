import sys
import subprocess
import os
import shutil
import numpy
type_to_label = {"SNV": 0, "Deletion": 1, "Insertion": 2, "Complex": 3, "MNV":3, "NONE": 4}

def format_data_item(jri, fisher):
    data = list()
    key = jri[2] + ":" + jri[3] # key is chrom:pos like "chr1:131022"
    #data.append(jri[3]) #start position
    #data.append(jri[4]) #end position   
    #jri[5] == cri[5] #refallele      
    #jri[6] == cri[6] #varallele      
    data.append(int(jri[7])) #totalposcoverage
    data.append(int(jri[8])) #positioncoverage
    data.append(int(jri[9])) #refForwardcoverage
    data.append(int(jri[10])) #refReversecoverage
    data.append(int(jri[11])) #varsCountOnForward
    data.append(int(jri[12])) #VarsCountOnReverse

    #jri[13] == cri[13] #genotype

    data.append(float(jri[14])) #frequency

    #jri[15] == cri[15] #strandbiasflag

    data.append(float(jri[16])) #meanPosition
    data.append(jri[17]) #pstd
    data.append(float(jri[18])) #meanQuality 
    data.append(jri[19]) #qstd
    index = 20
    if fisher:
        jri[index]  #pvalue
        index += 1
        jri[index]  #ratio
        index += 1
    
    data.append(float(jri[index])) #mapq
    index += 1
    data.append(float(jri[index])) #qratio
    index += 1
    data.append(float(jri[index])) #higreq
    index += 1
    data.append(float(jri[index])) #extrafreq
    index += 1
    data.append(int(jri[index])) #shift3
    index += 1
    data.append(float(jri[index])) #msi
    index += 1
    data.append(int(jri[index])) #msint
    index += 1
    data.append(float(jri[index])) #nm
    index += 1
    data.append(int(jri[index])) #hicnt
    index += 1
    data.append(int(jri[index])) #hicov
    index += 1

    #jri[30] == cri[30] #leftSequence
    #jri[31] == cri[31] #rightSequence
    #jri[32] == cri[32] #region
    #jri[33] == cri[33] #varType
    #jri[34]            # duprate
    if fisher:
        data.append(type_to_label[jri[35]])
    else:
        data.append(type_to_label[jri[33]])

    for i in range(len(data)):
        data[i] = float(data[i])

    return key, data

def read_strelka_data():
    data = list()
    key = "" + ":" + ""
    return key, data

def run_tools_and_get_data(fastvc_cmd, gen_cmd, strelka_cmd, base_path):
    tmpspace = os.path.join(base_path, "tmpspace") 
    if not os.path.exists(tmpspace):
        os.mkdir(tmpspace)
        os.mkdir(os.path.join(tmpspace, "strelka_space"))
        os.mkdir(os.path.join(tmpspace, "fastvc_space"))
    else:
        print("tmpspace exists! delete it first!")
        exit(-1)
    ret = subprocess.check_call(fastvc_cmd, shell = True)
    if not ret:
        print("fastvc runing error!!")
        exit(0)
    #--- read fastvc result file and format ---#
    fastvc_dict = dict()
    with open(os.path.join(base_path, "tmpspace/fastvc_space/out.txt"), 'r') as f:
        for line in f:
            items = line.split("\t")
            if len(items) == 36 :
                k, d = format_data_item(items, False)
                fastvc_dict[k] = d
            elif len(items) == 38 :
                k, d = format_data_item(items, True)
                fastvc_dict[k] = d
    #fastvc_data = numpy.asarray(fastvc_data)

    #--- generate strelka2 workspace and run ---#
    ret = subprocess.check_call(gen_cmd, shell = True)
    if ret:
        subprocess.check_call(strelka_cmd, shell = True)
    else:
        print("strelka gene workspace error!")
        exit(0)
    
    #--- read strelka2 result and format ---#
    tmp_res_sk2 = os.path.join(base_path, "tmpspace/strelka_space/results/results/variants/variants.vcf.gz")
    sk2_dict = list()
    with open(tmp_res_sk2, 'r') as f:
        for line in f:
            k, d = read_strelka_data(line)
            sk2_dict[k] = d 
    #sk2_data = numpy.asarray(sk2_data)
    #--- combine fastvc and sk2 result : all data merged into fastvc_dict---#
    fastvc_empty = [0.0 for i in range(23)]
    sk2_empty = [0.0 for i in range(14)]
    for k, v in fastvc_dict.items:
        if k not in sk2_dict:
           fastvc_dict[k] += sk2_empty
    for k, v in sk2_dict.items():
        if k in fastvc_dict:
            fastvc_dict[k] += v
        else:
            fastvc_dict[k] = fastvc_empty + v

    return fastvc_dict           


def prepare_cmds(fasta_file, region_file, bam_file, thread_number, base_path):
    #--- fastvc cmd prepareing ---#
    fvc_list = list() 
    fastvc_path = ""
    fvc_list.append(fastvc_path)
    fvc_list.append("-i {}".format(region_file))
    fvc_list.append("-G {}".format(fasta_file))
    fvc_list.append("-f 0.01")
    fvc_list.append("-N NA12878")
    fvc_list.append("-b {}".format(bam_file))
    fvc_list.append("-c 1 -S 2 -E 3 -g 4")
    fvc_list.append("--fisher")
    fvc_list.append("--th {}".format(thread_number))
    fvc_list.append("--out {}".format(os.path.join(base_path, "tmpspace/fastvc_space/out.txt")))
    fastvc_cmd = " ".join(fvc_list)    

    #--- strelka cmd prepareing ---#
    #-- 1. generate workspace and script --#
    sk2_conf_path = ""
    gen_cmd = "{} --bam {} --referenceFasta {} --callRegions {} --runDir {}".format(sk2_conf_path, 
        bam_file, fasta_file, region_file, os.path.join(base_path, "tmpspace/strelka_space")) 

    #-- 2. strelka run command --#
    sk2_cmd = "{}/tmpspace/strelka_space/runWorkflow.py  -m local -j {}".format(base_path, thread_number)

    return fastvc_cmd, gen_cmd, sk2_cmd

if __name__ == "__main__":
    region_file = "/home/old_home/haoz/workspace/data/NA12878/ConfidentRegions.bed"
    fasta_file = "/home/old_home/haoz/workspace/data/hg38/hg38.fa"
    bam_file = "/home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam"
    base_path = "/home/old_home/haoz/workspace"
    fastvc_cmd, gen_cmd, sk2_cmd = prepare_cmds(fasta_file, region_file, bam_file, 40, base_path)
    #merged_data: dict: key=chrom:pos, value = [fastvc_feature, sk2_feature] (23 + 14 dim)
    merged_data = run_tools_and_get_data(fastvc_cmd, gen_cmd, sk2_cmd, base_path) 