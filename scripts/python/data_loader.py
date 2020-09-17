import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import torch
from pandas import DataFrame
from sklearn import preprocessing 
from sklearn.model_selection import train_test_split
import pickle
from features import features_to_index as fe2i, fvc_selected_features as fvc_sf

type_to_label = {"SNV":     [1, 0, 0, 0, 0, 0], 
                "Deletion": [0, 1, 0, 0, 0, 0], 
                "Insertion":[0, 0, 1, 0, 0, 0], 
                "Complex":  [0, 0, 0, 1, 0, 0], 
                "MNV":      [0, 0, 0, 0, 1, 0], 
                "NONE":     [0, 0, 0, 0, 0, 1]}


FVC_FEATURES = 19
SK2_FEATURES = 13
#SK2_FEATURES = 0

def origional_format_data_item(jri, fisher):
    data = list()
    # key is chrom:pos like "chr1:131022:A:AT"
    key = jri[2] + ":" + jri[3] + ":" + jri[5] + ":" + jri[6]
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
        data.extend(type_to_label[jri[35]])
    else:
        data.extend(type_to_label[jri[33]])

    if len(data) != FVC_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        #print("origin data:" , jri)
    return key, data

def format_data_item(jri, fisher = True):
    #fisher should always be true
    data = list()
    # key is chrom:pos like "chr1:131022:A:AT"
    key = jri[2] + ":" + jri[3] + ":" + jri[5] + ":" + jri[6]
    data.append(len(jri[fe2i["RefAllel"]])) #refallele len
    data.append(len(jri[fe2i["VarAllel"]])) #varallel len
    for sf in fvc_sf:
        data.append(jri[fe2i[sf]])
    
    if fisher:
        data.extend(type_to_label[jri[fe2i["varType"]]])
    else:
        data.extend(type_to_label[jri[fe2i["varType"]]])

    if len(data) != FVC_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        #print("origin data:" , jri)
    return key, data

def read_strelka_data(items):
    data = list()
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
    #GT  : GQ : GQX : DP : DPF : AD  : ADF : ADR : SB  : FT              : PL
    #0/1 : 3  : 0   : 1  : 0   : 0,1 : 0,0 : 0,1 : 0.0 : LowGQX;LowDepth : 28,3,0

    #CIGAR=1M2D     ;RU=TG   ;REFREP=3    ;IDREP=2    ;MQ=35 
    #CIGAR=1M2I,1M4I;RU=GT,GT;REFREP=12,12;IDREP=13,14;MQ=59
    #GT  : GQ : GQX : DPI : AD     : ADF  : ADR  : FT   :PL 
    #0/1 : 35 : 7   : 15  : 12,3   : 10,2 : 2,1  : PASS :32,0,240
    #1/2 : 196: 19  : 44  : 0,18,14: 0,8,9:0,10,5: PASS :671,266,199,360,0,301
    #key = items[0] + ":" + items[1]
    #return key, [0] * SK2_FEATURES
    ALTs = items[4].split(',')
    key = items[0] + ":" + items[1] + ":" + items[3] + ":" + ALTs[0]
    mult_alt = False
    if len(ALTs) > 1:
        mult_alt = True
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
    #data.append(0 if info_dict['GT'] == '.'  else 1 + int(info_dict['GT'][0]) + int(info_dict['GT'][-1]))
    data.append(0 if info_dict['GQ'] == '.' else info_dict['GQ']) #GQ
    data.append(0 if info_dict['GQX'] == '.' else info_dict['GQX'])#GQX
    #DP(SNV) or DPI(indels)
    if is_indel:
        data.append(info_dict['DPI'])
    else:
        data.append(info_dict['DP'])
    #注：因为可能出现两个ALT的情况，这里为了测试暂时只取第一个ALT
    data.extend([0, 0] if info_dict['AD'] == '.' else info_dict['AD'].split(',')[:2]) #AD 2
    #data.extend([0, 0] if info_dict['ADF'] == '.' else info_dict['ADF'].split(',')[:2]) #ADF 2
    #data.extend([0, 0] if info_dict['ADR'] == '.' else info_dict['ADR'].split(',')[:2]) #ADR 2
    #filter: one-hot type, [1,0] if pass, else [0,1]
    data.extend([1, 0]if info_dict['FT'] == "PASS" else [0, 1])
    #PL 3
    if len(info_dict['GT']) == 1:
        data.extend([0, 0] if info_dict['PL'] == '.' else info_dict['PL'].split(',')[:2])
        data.extend([0])
    else:
        data.extend([0, 0, 0] if info_dict['PL'] == '.' else info_dict['PL'].split(',')[:3])

    if len(data) != SK2_FEATURES:
        print("sk2 data length error: \n", len(data), data, "ori: \n", items)
    return key, data

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
    #'''------------------------先不要strelka了---------------------------------------
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
    for k, v in fastvc_dict.items():
        if k not in sk2_dict:
           fastvc_dict[k] += sk2_empty
    for k, v in sk2_dict.items():
        if k in fastvc_dict:
            fastvc_dict[k] += v
        else:
            fastvc_dict[k] = fastvc_empty + v
   #--------------------------------------------------------------- ''' 
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
                    site = chrom + ":" + pos + ":" + ref + ":" + alt.split(",")[0]
                    truth_vars.add(site)
                    #truth_vars[site] = list([ref, alt])  
    labels_dict = {}
    positive_num = 0
    negtive_num = 0
    for k, v in data_dict.items():
        if k in truth_vars:
            #labels_dict[k] = [1, 0]
            labels_dict[k] = 1
            positive_num += 1
        else:
            #labels_dict[k] = [0, 1]
            labels_dict[k] = 0
            negtive_num += 1
    return positive_num, negtive_num, labels_dict

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
    def __init__(self, data):
        self.inputs = data[0]
        self.labels = data[1]

    def __getitem__(self, index):
        input, label = self.inputs[index], self.labels[index]
        return input, np.asarray(label)

    def __len__(self):
        return len(self.labels)

class Dataset:
    def __init__(self, reload, re_exec, pama_list, base_path, truth_path):
        self.data = dict()
        if not reload:
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
            pos_num, neg_num, merged_label_dict = get_labels_dict(merged_data_dict, truth_path)
            print("get label done, size: {}, pos_num: {}, neg_num: {}".format(len(merged_label_dict), pos_num, neg_num))
            keys = list()
            self.inputs = list()
            self.labels = list()
            for k, v in merged_data_dict.items():
                #self.data.append([k, np.asarray(v), merged_label_dict[k]])
                keys.append(k)
                self.inputs.append(v)
                self.labels.append(merged_label_dict[k])
            #--- standlizaton ---#            
            '''
            min_max_scaler = preprocessing.MinMaxScaler()
            self.inputs = min_max_scaler.fit_transform(self.inputs)
            '''
            #---inputs Normalization ---#
            self.inputs = np.asfarray(self.inputs)
            print("start normalization...")
            self.inputs = preprocessing.normalize(self.inputs, axis = 0, norm = 'max') 
            '''
            np.set_printoptions(precision = 3, threshold=1000)
            for i in range(10000):
                print(self.keys[i], self.inputs[i], self.labels[i])
            exit(0)
            '''
            print("normalization done")
            print("FastvcDataset init over")
            #self.df = DataFrame([self.keys, self.inputs, self.labels], columns=["keys", "inputs", "labels"])
    def split(self, test_size = 0.2, random_state = 0):
        print("spliting data...")
        x_train, x_test, y_train, y_test = train_test_split(self.inputs, self.labels, 
                                test_size = test_size, random_state = random_state)

        #x_train = x_train[0:600000] 
        #y_train = y_train[0:600000] 
        #x_test = x_test[0:60000] 
        #y_test = y_test[0:60000] 
        print("train summary: toatl:{}, 0:{}, 1:{}".format(len(y_train), len(y_train) - sum(y_train), sum(y_train)) )
        print("test summary: toatl:{}, 0:{}, 1:{}".format(len(y_test), len(y_test) - sum(y_test), sum(y_test)) )
        self.data['train'] = [x_train, y_train]
        self.data['test'] = [x_test, y_test]

    def load(self, store_path):
        if not os.path.exists(store_path):
            print("file {} do not exists!".format(store_path))
            exit(-1)
        with open(store_path, 'rb') as f:
            self.data = pickle.load(f)
        print("load data over!")
        
    def store(self, store_path):
        if os.path.exists(store_path):
            print("file {} exists!".format(store_path))
            exit(-1)
        with open(store_path, 'wb') as f:
            pickle.dump(self.data, f)
        print("dump finished!")