import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import torch
from pandas import DataFrame
from sklearn.preprocessing import normalize

type_to_label = {"SNV": 0, "Deletion": 1, "Insertion": 2, "Complex": 3, "MNV":3, "NONE": 4}

FVC_FEATURES = 26
SK2_FEATURES = 0 

def format_data_item(jri, fisher):
    data = list()
    key = jri[2] + ":" + jri[3] # key is chrom:pos like "chr1:131022"
    #data.append(jri[3]) #start position
    #data.append(jri[4]) #end position   
    #jri[5] == cri[5] #refallele      
    #jri[6] == cri[6] #varallele      
    data.append(len(jri[5])) #refallele len
    data.append(len(jri[6])) #varallel len
    data.append(jri[7]) #totalposcoverage
    data.append(jri[8]) #positioncoverage
    data.append(jri[9]) #refForwardcoverage
    data.append(jri[10]) #refReversecoverage
    data.append(jri[11]) #varsCountOnForward
    data.append(jri[12]) #VarsCountOnReverse
    #jri[13] == cri[13] #genotype
    data.append(jri[14]) #frequency
    #jri[15] == cri[15] #strandbiasflag
    data.append(jri[16]) #meanPosition
    data.append(jri[17]) #pstd
    data.append(jri[18]) #meanQuality 
    data.append(jri[19]) #qstd
    index = 20
    if fisher:
        data.append(jri[index])  #pvalue
        index += 1
        data.append(jri[index])  #ratio
        index += 1
    
    data.append(jri[index]) #mapq
    index += 1
    data.append(jri[index]) #qratio
    index += 1
    data.append(jri[index]) #higreq
    index += 1
    data.append(jri[index]) #extrafreq
    index += 1
    data.append(jri[index]) #shift3
    index += 1
    data.append(jri[index]) #msi
    index += 1
    data.append(jri[index]) #msint
    index += 1
    data.append(jri[index]) #nm
    index += 1
    data.append(jri[index]) #hicnt
    index += 1
    data.append(jri[index]) #hicov
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

def read_strelka_data(items):
    data = list()
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
    #GT  : GQ : GQX : DP : DPF : AD  : ADF : ADR : SB  : FT              : PL
    #0/1 : 3  : 0   : 1  : 0   : 0,1 : 0,0 : 0,1 : 0.0 : LowGQX;LowDepth : 28,3,0

    #CIGAR=1M2D;RU=TG;REFREP=3;IDREP=2;MQ=35 
    #GT  : GQ : GQX : DPI : AD   : ADF  : ADR  : FT   :PL 
    #0/1 : 35 : 7   : 15  : 12,3 : 10,2 : 2,1  : PASS :32,0,240
    key = items[0] + ":" + items[1]
    return key, [0] * SK2_FEATURES
    '''
    #print(items)
    key = items[0] + ":" + items[1]
    data.append(0 if items[5] == '.' else items[5]) #QUAL

    INFO = items[7]
    FORMAT = items[8].split(':')
    VALUES = items[9].split(':')
    info_dict = {}
    for i in range(len(FORMAT)):
        info_dict[FORMAT[i].strip()] = VALUES[i].strip()
    MQ = INFO.split(";")[-1] 
    if MQ[0:3] == "MQ=":
        data.append(MQ[3:]) #map quality
    else:
        print("error! invalide data(MQ)!", MQ)
        exit(-1)
    is_indel = 0
    if INFO.split(":")[0][:5] == "CIGAR":
        is_indel = 1
    data.append(is_indel) #if this var is indel
    #genotype 1: 0/0, 2: 0/1, 3: 1/1, 0: .
    data.append(0 if info_dict['GT'] == '.'  else 1 + int(info_dict['GT'][0]) + int(info_dict['GT'][-1]))
    data.append(0 if info_dict['GQ'] == '.' else info_dict['GQ']) #GQ
    data.append(0 if info_dict['GQX'] == '.' else info_dict['GQX'])#GQX
    #DP(SNV) or DPI(indels)
    if is_indel:
        data.append(info_dict['DPI'])
    else:
        data.append(info_dict['DP'])
    data.extend([0, 0] if info_dict['AD'] == '.' else info_dict['AD'].split(',')) #AD 2
    data.extend([0, 0] if info_dict['ADF'] == '.' else info_dict['ADF'].split(',')) #ADF 2
    data.extend([0, 0] if info_dict['ADR'] == '.' else info_dict['ADR'].split(',')) #ADR 2
    #filter
    data.append(1 if info_dict['FT'] == "PASS" else 0)
    #PL 3
    data.extend([0, 0, 0] if info_dict['PL'] == '.' else info_dict['PL'].split(','))
    #print(len(data),data)
    if len(data) != SK2_FEATURES:
        #print(items)
        #print("data length error", len(data), data)
        data = [0] * SK2_FEATURES
    #else:
    #    print("normaldata")
    try:
        data = list(map(float, data))
    except:
        print("error data:", data)
    return key, data
    '''

def get_data(fvc_result_path, sk2_result_path):
    #--- read fastvc result file and format ---#
    fastvc_dict = dict()
    with open(fvc_result_path, 'r') as f:
        for line in f:
            items = line.split("\t")
            if len(items) == 36:
                k, d = format_data_item(items, False)
                fastvc_dict[k] = d
            elif len(items) == 38 :
                k, d = format_data_item(items, True)
                fastvc_dict[k] = d
    print("get fastvc data done: ", len(fastvc_dict))
    #--- read strelka2 result and format ---#
    sk2_dict = dict()
    with open(sk2_result_path, 'r') as f:
        for line in f:
            if line[0] == '#': 
                continue
            items = line.split('\t')
            if len(items) == 10:
                k, d = read_strelka_data(items)
                sk2_dict[k] = d 
    print("get sk2 data done: ", len(sk2_dict))
    #--- combine fastvc and sk2 result : all data merged into fastvc_dict---#
    fastvc_empty = [0.0] * FVC_FEATURES
    sk2_empty = [0.0] * SK2_FEATURES
    lens = {24:0, 26:0, 16:0, 18:0}
    '''
    for k, v in fastvc_dict.items():
        lens[len(v)] += 1
        if k not in sk2_dict:
           fastvc_dict[k] += sk2_empty
    for k, v in sk2_dict.items():
        lens[len(v)] += 1
        if k in fastvc_dict:
            fastvc_dict[k] += v
        else:
            fastvc_dict[k] = fastvc_empty + v
    '''
    print("fvc lens: ", lens)
    return fastvc_dict           

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
        exit(-1)
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
    sk2_relative_varpath = "tmpspace/strelka_space/results/results/variants/variants.vcf.gz"
    tmp_res_sk2 = os.path.join(base_path, sk2_relative_varpath)
    sk2_dict = list()
    with open(tmp_res_sk2, 'r') as f:
        for line in f:
            k, d = read_strelka_data(line)
            sk2_dict[k] = d 
    #sk2_data = numpy.asarray(sk2_data)
    #--- combine fastvc and sk2 result : all data merged into fastvc_dict---#
    fastvc_empty = [0.0 for i in range(FVC_FEATURES)]
    sk2_empty = [0.0 for i in range(SK2_FEATURES)]
    for k, v in fastvc_dict.items:
        if k not in sk2_dict:
           fastvc_dict[k] += sk2_empty
    for k, v in sk2_dict.items():
        if k in fastvc_dict:
            fastvc_dict[k] += v
        else:
            fastvc_dict[k] = fastvc_empty + v

    return fastvc_dict           

def get_labels_dict(data_dict, truth_path):
    #truth_vars = dict()
    truth_vars = set()
    with open(truth_path, 'r') as f:
        for var in f:
            items = var.split('\t')
            if(len(items) == 10 ):
                chrom, pos, id, ref, alt, _, filter = items[:7]         
                #if len(chrom) < 6 and filter == "PASS" and (len(ref) > 1 or len(alt) > 1) :
                if len(chrom) < 6 and filter == "PASS":
                    site = chrom + ":" + pos
                    truth_vars.add(site)
                    #truth_vars[site] = list([ref, alt])  
    labels_dict = {}
    for k, v in data_dict.items():
        if k in truth_vars:
            labels_dict[k] = [1, 0]
        else:
            labels_dict[k] = [0, 1]
    return labels_dict

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

class FastvcDataset(torch.utils.data.Dataset):
    #def __init__(self, region_file, fasta_file, bam_file, base_path, truth_path):
    def __init__(self, re_exec, pama_list, base_path, truth_path):
        merged_data_dict = {}
        if re_exec:
            region_file, fasta_file, bam_file = pama_list
            fastvc_cmd, gen_cmd, sk2_cmd = prepare_cmds(fasta_file, region_file, bam_file, 40, base_path)
            #merged_data: dict: key=chrom:pos, value = [fastvc_feature, sk2_feature] (FVC_FEATURE2 + SK2_FEATURES dim)
            merged_data_dict = run_tools_and_get_data(fastvc_cmd, gen_cmd, sk2_cmd, base_path) 
        else:
            fvc_res_path, sk2_res_path = pama_list
            merged_data_dict = get_data(fvc_res_path, sk2_res_path) 
        assert(len(merged_data_dict) > 0)

        print("get merged data done, merged data dict size: ", len(merged_data_dict))
        merged_label_dict = get_labels_dict(merged_data_dict, truth_path)
        print("get label done, size:", len(merged_label_dict))
        self.keys = list()
        self.inputs = list()
        self.labels = list()
        for k, v in merged_data_dict.items():
            #self.data.append([k, np.asarray(v), merged_label_dict[k]])
            self.keys.append(k)
            self.inputs.append(v)
            self.labels.append(merged_label_dict[k])
        #---inputs Normalization ---#
        self.inputs = np.asfarray(self.inputs)
        print("start normalization...")
        self.inputs = normalize(self.inputs, axis = 0, norm = 'l2') 
        print("FastvcDataset init over")

    def __getitem__(self, index):
        key, input, label = self.keys[index], self.inputs[index], self.labels[index]
        return key, input, np.asarray(label)

    def __len__(self):
        return len(self.keys)
