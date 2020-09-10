import torch
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torchvision import transforms
import torchvision
import datetime
import os

from nn_net import Net
from data_loader import FastvcDataset
import data_loader 

def train_epoch(train_loader, net, optimizer):
    epoch_loss = [] 
    runing_loss = 0.0
    use_cuda = False
    for i, data in enumerate(train_loader, 0):
        (key, inputs, labels) = data #TODO #DONE
        inputs, labels = Variable(inputs).float(), Variable(labels).float()
        if use_cuda:
            inputs, labels = inputs.cuda(), labels.cuda()
        optimizer.zero_grad()
        outputs = net(inputs) #outputs is the prob of each class(P or N)
        loss_ = loss_func(outputs, labels)
        epoch_loss.append(loss_.data)
        runing_loss += loss_.data
        loss_.backward()

        #update wight
        optimizer.step()

        if i % 2000 == 1999:
            print("[%5d] loss: %.3f" % (i + 1, loss_.item()))
            runing_loss = 0.0

    return epoch_loss


#--------------------------------------------------------#
region_file = "/home/old_home/haoz/workspace/data/NA12878/ConfidentRegions.bed"
fasta_file = "/home/old_home/haoz/workspace/data/hg38/hg38.fa"
bam_file = "/home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam"
base_path = "/home/old_home/haoz/workspace/FastVC/workspace"
truth_path =  "/home/old_home/haoz/workspace/data/NA12878/vcfs/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf"
out_dir = os.path.join(base_path, "out")
re_exec = False
strelka2_result_path = "/home/old_home/haoz/workspace/VCTools/strelka-2.9.10.centos6_x86_64/hg38run_germline_test/results/variants/variants.vcf"
fastvc_result_path = "/home/old_home/haoz/workspace/FastVC/detection_result/NA12878/of.txt"
#--------------------------------------------------------#
n_feature = data_loader.FVC_FEATURES + data_loader.SK2_FEATURES # fastvc 25 + sk2 14
net = Net(n_feature, 50, 2)
max_epoch = 100
save_freq = 20 # save every xx save_freq
n_epoch = 0
batch_size = 64
nthreads = 20
optimizer = optim.SGD(net.parameters(), lr = 0.03)

loss_func = torch.nn.MSELoss()

if re_exec:
    train_set = FastvcDataset(re_exec, [region_file, fasta_file, bam_file], 
                            base_path, truth_path)
else:
    train_set = FastvcDataset(re_exec, [fastvc_result_path, strelka2_result_path],
                            base_path, truth_path)
train_loader = torch.utils.data.DataLoader(train_set, 
                            batch_size = batch_size, 
                            num_workers = nthreads, pin_memory = True)

for epoch in range(max_epoch):
    print("epoch", epoch, " processing....")
    n_epoch += 1
    epoch_loss = train_epoch(train_loader, net, optimizer)
    print("mean loss of epoch %d is: %f" % (epoch, sum(epoch_loss) / len(epoch_loss)))
    if n_epoch % save_freq == 0:
        tag = "fastvc_{}".format(datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
        torch.save({
            "state_dict": net.state_dict(),
            "tag": tag,
            "epoch": n_epoch,
            }, '{}/models/checkpoint_{}_ecpch{}.pth'.format(out_dir, tag, n_epoch))

tag = "fastvc_{}".format(datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
torch.save({
    "state_dict": net.state_dict(),
    "tag": tag,
    "epoch": n_epoch,
    }, '{}/models/checkpoint_{}_ecpch{}.pth'.format(out_dir, tag, n_epoch))

print("training done!")
print("final model:", '{}/models/checkpoint_{}_ecpch{}.pth'.format(out_dir, tag, n_epoch))
