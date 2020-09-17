import torch 
from torch.autograd import Variable 
import torch.nn as nn
import torch.nn.functional as F 
import matplotlib.pyplot as plt 
  
  
class Net(torch.nn.Module): 
  def __init__(self, n_features, hidden_sizes, n_output): 
    super(Net, self).__init__()  
    '''
    self.layers = []
    self.layers.append(nn.Linear(n_features, hidden_sizes[0]))
    for i in range(len(hidden_sizes) - 1):
      print("add hidden layer:({} x {})".format(hidden_sizes[i], hidden_sizes[i+1]))
      self.layers.append(nn.Linear(hidden_sizes[i], hidden_sizes[i + 1]))
    '''
    self.hidden0 = nn.Linear(n_features, hidden_sizes[0])
    assert(len(hidden_sizes) == 4)
    self.hidden1 = nn.Linear(hidden_sizes[0], hidden_sizes[1])
    self.hidden2 = nn.Linear(hidden_sizes[1], hidden_sizes[2])
    self.hidden3 = nn.Linear(hidden_sizes[2], hidden_sizes[3])
    self.predict = nn.Linear(hidden_sizes[-1], n_output) # output layer
  
  def forward(self, x): 
    '''
    for i in range(len(self.layers)):
      x = F.relu(self.layers[i](x))
    '''
    #print("----------------print out input--------------")
    #print("1\n", x)
    x = torch.relu(self.hidden0(x)) 
    #print("2\n",x)
    x = torch.relu(self.hidden1(x)) 
    #print("3\n",x)
    x = torch.relu(self.hidden2(x))
    #print("4\n",x)
    x = torch.relu(self.hidden3(x))
    #print("5\n",x)
    x = self.predict(x) 
    #print("6\n", x)
    return x
    #return F.softmax(x, dim = 1)

  def initialize_weights(self):
    print("initializeing net work")
    for m in self.modules():
      print(m)
      if isinstance(m, nn.Linear):
        # print(m.weight.data.type())
        # input()
        # m.weight.data.fill_(1.0)
        nn.init.xavier_uniform_(m.weight, gain=nn.init.calculate_gain('relu'))
        #nn.init.xavier_uniform_(m.weight)

'''  
opitmizer = torch.optim.SGD(net.parameters(),lr=0.03)
loss_fun = nn.MSELoss()   #选择 均方差为误差函数

# 定义优化器和损失函数 
optimizer = torch.optim.SGD(net.parameters(), lr=0.05) # 传入网络参数和学习率 
loss_function = torch.nn.MSELoss() # 最小均方误差 
  
# 神经网络训练过程 
plt.ion()  # 动态学习过程展示 
plt.show() 
  
for t in range(300): 
  prediction = net(x) # 把数据x喂给net，输出预测值 
  loss = loss_function(prediction, y) # 计算两者的误差，要注意两个参数的顺序 
  optimizer.zero_grad() # 清空上一步的更新参数值 
  loss.backward() # 误差反相传播，计算新的更新参数值 
  optimizer.step() # 将计算得到的更新值赋给net.parameters() 
  
  # 可视化训练过程 
  if (t+1) % 10 == 0: 
    plt.cla() 
    plt.scatter(x.data.numpy(), y.data.numpy()) 
    plt.plot(x.data.numpy(), prediction.data.numpy(), 'r-', lw=5) 
    plt.text(0.5, 0, 'L=%.4f' % loss.data[0], fontdict={'size': 20, 'color': 'red'}) 
    plt.pause(0.1)
'''