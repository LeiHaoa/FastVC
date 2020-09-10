import torch 
from torch.autograd import Variable 
import torch.nn.functional as F 
import matplotlib.pyplot as plt 
  
  
class Net(torch.nn.Module): 
  def __init__(self, n_feature, n_hidden, n_output): 
    super(Net, self).__init__()  
    self.hidden1 = torch.nn.Linear(n_feature, n_hidden) # hidden layer 
    #self.hidden2 = torch.nn.Linear(n_hidden, n_hidden) # hidden layer 
    #dropout???
    self.predict = torch.nn.Linear(n_hidden, n_output) # output layer
  
  def forward(self, x): 
    x = F.relu(self.hidden1(x)) 
    #x = F.relu(self.hidden2(x))
    x = self.predict(x) 
    return F.softmax(x, dim = 1)

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