import datetime
import os
import sys

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torchvision
from torch.autograd import Variable
from torchvision import transforms

import data_loader
from data_loader import Dataset, FastvcDataset
from MyLogger import Logger
from nn_net import Net

sys.stdout = Logger(filename = "./mcall.out")
#--------------------------------------------------------#
region_file = "/home/old_home/haoz/workspace/data/NA12878/ConfidentRegions.bed"
fasta_file = "/home/old_home/haoz/workspace/data/hg38/hg38.fa"
bam_file = "/home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam"
base_path = "/home/old_home/haoz/workspace/FastVC/workspace"
truth_path =  "/home/old_home/haoz/workspace/data/NA12878/vcfs/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf"
out_dir = os.path.join(base_path, "out/models")
re_exec = False
strelka2_result_path = "/home/old_home/haoz/workspace/VCTools/strelka-2.9.10.centos6_x86_64/hg38run_40th/results/variants/variants.vcf"
#strelka2_result_path = "/home/old_home/haoz/workspace/VCTools/strelka-2.9.10.centos6_x86_64/hg38run_germline_test/results/variants/variants.vcf"
fastvc_result_path = "/home/old_home/haoz/workspace/FastVC/detection_result/NA12878/out_fisher.vcf"
#--------------------------------------------------------#

reload_from_dupfile = False #load from file(True) or compute generate data again(Fasle)
data_path = "./call_dataset.pkl"
if re_exec:
    dataset = Dataset(reload_from_dupfile, re_exec, [region_file, fasta_file, bam_file], 
                            base_path, truth_path)
else:
    dataset = Dataset(reload_from_dupfile, re_exec, [fastvc_result_path, strelka2_result_path],
                            base_path, truth_path)
if reload_from_dupfile:
    dataset.load(data_path)
else:
    if os.path.exists(data_path):
        os.remove(data_path)
    dataset.split(random_state = None)
    dataset.store(data_path)
#------------------------network setting---------------------#
n_feature = data_loader.FVC_FEATURES + data_loader.SK2_FEATURES # fastvc 31 + sk2 17
net = Net(n_feature, [40, 60, 60, 10] , 2)
use_cuda = False
batch_size = 32 
nthreads = 20
#------------------------------------------------------------#
 
epoch_loss = [] 
runing_loss = 0.0
use_cuda = False
#optimizer.zero_grad()
total = 0
truth = 0
false = 0
tmp_file = open("./tmp.out", 'w')
for i, data in enumerate(test_loader, 0):
    inputs, labels = data #TODO #DONE
    inputs = Variable(inputs).float()
    total += len(inputs)
    if use_cuda:
        inputs = inputs.cuda()

    outputs = net(inputs) #outputs is the prob of each class(P or N)
    _, predicted = torch.max(outputs, 1)
    compare_labels = (predicted == labels)
    false_preds = np.where(compare_labels.numpy() == 0)[0]
    false += len(false_preds)
print("test epoch: [total]:{}, [false]:{}, [truth]:{}, error rate:{}".format(total, false, total - false, false/total) )
tmp_file.close()