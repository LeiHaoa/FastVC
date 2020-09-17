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

sys.stdout = Logger()
def train_epoch(train_loader, net, optimizer):
    epoch_loss = [] 
    runing_loss = 0.0
    use_cuda = False
    #optimizer.zero_grad()
    for i, data in enumerate(train_loader, 0):
        inputs, labels = data #TODO #DONE
        inputs, labels = Variable(inputs).float(), Variable(labels).long()
        if use_cuda:
            inputs, labels = inputs.cuda(), labels.cuda()
        optimizer.zero_grad()
        outputs = net(inputs) #outputs is the prob of each class(P or N)
        loss_ = loss_func(outputs, labels)
        loss_.backward()
        #update wight
        optimizer.step()

        epoch_loss.append(loss_.item())
        runing_loss += loss_.item()

        if i % 10000 == 9999:
            print("[%5d] loss: %.3f" % (i + 1, runing_loss / 10000))
            #np.set_printoptions(precision = 3, threshold=10000)
            #print("[inputs info]: \n {}".format(inputs.data.cpu().numpy()))
            #print("[pred info]:\n {}".format(outputs.data.cpu().numpy()))
            #print("[label info]:\n {}".format(labels.data.cpu().numpy()))
            #print('[net info]:')
            #for layer in net.modules():
            #    if isinstance(layer, nn.Linear):
            #        print(layer.weight)
            runing_loss = 0.0

    return epoch_loss

def write_to_tmp_file(p, l, f):
    for i in range(len(p)):
        f.write(str(p.numpy()[i]) + " - " + str(l.numpy()[i]) + '\n')

def test_epoch(test_loader, net, optimizer):
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
        write_to_tmp_file(predicted, labels, tmp_file) #----debug
        compare_labels = (predicted == labels)
        false_preds = np.where(compare_labels.numpy() == 0)[0]
        false += len(false_preds)
    print("test epoch: [total]:{}, [false]:{}, [truth]:{}, error rate:{}".format(total, false, total - false, false/total) )
    tmp_file.close()

    return epoch_loss

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

reload_from_dupfile = True #load from file(True) or compute generate data again(Fasle)
data_path = "./dataset.pkl"
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
max_epoch = 100 
use_cuda = False
save_freq = 20 # save every xx save_freq
n_epoch = 0
batch_size = 32 
nthreads = 20
#init
#net.initialize_weights()
#optimizer
#optimizer = optim.SGD(net.parameters(), lr = 0.1, momentum = 0.9)
optimizer = optim.Adam(net.parameters(), lr = 0.01)
#loss_func = torch.nn.MSELoss()
loss_func = torch.nn.CrossEntropyLoss()
#loss_func = torch.nn.BCELoss()
#------------------------------------------------------------#

for epoch in range(max_epoch):
    print("epoch", epoch, " processing....")
    if (not reload_from_dupfile) and (epoch != 0):
        dataset.split(random_state = None)
    test_dataset = FastvcDataset(dataset.data['test'])
    test_loader = torch.utils.data.DataLoader(test_dataset, 
                                batch_size = batch_size, shuffle = True,
                                num_workers = nthreads, 
                                pin_memory = True)
    train_dataset = FastvcDataset(dataset.data['train'])
    train_loader = torch.utils.data.DataLoader(train_dataset, 
                                batch_size = batch_size, shuffle = True,
                                num_workers = nthreads, 
                                pin_memory = True)

    n_epoch += 1
    #epoch_loss = train_epoch(train_loader, net, optimizer)
    #------------------------train epoch-----------------
    epoch_loss = [] 
    runing_loss = 0.0
    #optimizer.zero_grad()
    for i, data in enumerate(train_loader, 0):
        inputs, labels = data #TODO #DONE
        inputs, labels = Variable(inputs).float(), Variable(labels).long()
        if use_cuda:
            inputs, labels = inputs.cuda(), labels.cuda()
        optimizer.zero_grad()
        outputs = net(inputs) #outputs is the prob of each class(P or N)
        loss_ = loss_func(outputs, labels)
        loss_.backward()
        #update wight
        optimizer.step()

        epoch_loss.append(loss_.item())
        runing_loss += loss_.item()

        if i % 10000 == 9999:
            print("[%5d] loss: %.3f" % (i + 1, runing_loss / 9999))
            #np.set_printoptions(precision = 3, threshold=10000)
            #print("[inputs info]: \n {}".format(inputs.data.cpu().numpy()))
            #print("[pred info]:\n {}".format(outputs.data.cpu().numpy()))
            #print("[label info]:\n {}".format(labels.data.cpu().numpy()))
            #print('[net info]:')
            #for layer in net.modules():
            #    if isinstance(layer, nn.Linear):
            #        print(layer.weight)
            runing_loss = 0.0

    print("mean loss of epoch %d is: %f" % (epoch, sum(epoch_loss) / len(epoch_loss)))
    test_epoch(test_loader, net, optimizer)
    if n_epoch % save_freq == (save_freq - 1):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        tag = "fastvc_{}".format(datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
        torch.save({
            "state_dict": net.state_dict(),
            "tag": tag,
            "epoch": n_epoch,
            }, '{}/checkpoint_{}_ecpch{}.pth'.format(out_dir, tag, n_epoch))

if not os.path.exists(out_dir):
    os.mkdir(out_dir)
tag = "fastvc_{}".format(datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
torch.save({
    "state_dict": net.state_dict(),
    "tag": tag,
    "epoch": n_epoch,
    }, '{}/checkpoint_{}_ecpch{}.pth'.format(out_dir, tag, n_epoch))
print("training done!")
print("final model:", '{}/models/checkpoint_{}_ecpch{}.pth'.format(out_dir, tag, n_epoch))
