import numpy as np
import os 
import subprocess

def read_data(file_name):
    return np.loadtxt(file_name)

data_list = ["./para_data/fin_data/" + name for name in os.listdir("./para_data/fin_data/")]
para_data = []
res_data = []
for data in data_list:
    this_res = read_data(data)
    para_data.append(this_res[0,:])
    res_data.append(this_res[1,:])
para_data = np.array(para_data)
res_data = np.array(res_data)
event_number = np.shape(res_data)[0]
count = 0
for i in range(0,event_number):
    if (i == 1000):
        break
    print("Doing iter in event {}".format(i))
    this_res = res_data[i,:]
    this_para = para_data[i,:]
    try_time = 0
    while(try_time < 20):
        try_time+=1
        print("this is try time {}".format(try_time))
        new_para = np.random.uniform(0.5,1.5,4)*this_para
        input_2 = ' '.join(str(x) for x in new_para)
        input_1 = ' '.join(str(x) for x in this_res)
        back = os.popen("./run.exe " + input_1 + " " + input_2)
        tmp_a = back.read()
        if(tmp_a == "error of nan in bessel"):
            print("nan in bessel, try again...")
            continue
        else:
            break
    if(try_time == 20 and (tmp_a == "error of nan in bessel")):
        print("failed, next event")
        continue
    print("real para  lambda_e , lambda_Pi, gamma1, gamma2 = {}".format(' '.join(str(x) for x in this_para)))
    print("In cal res lambda_e , lambda_Pi, gamma1, gamma2 = {}".format(tmp_a))
    array = np.fromstring(tmp_a, sep=' ')
    err_sum = np.max((abs(array - this_para)/abs(this_para+1e-7)))
    print('maximum relative error = {:.2f}%'.format(err_sum*100))
    print("\n")
    if err_sum < 1e-1:
        count+=1
print(count/1000)
