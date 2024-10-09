import os

dir_path = '/home/yuyu/bike/LDPC-BIKE/bf_data1_1k'

files = [f for f in os.listdir(dir_path) if f.startswith('output')]

count = {}

for filename in files:
    file_path = os.path.join(dir_path, filename)
    
    t = int(filename.split('_')[-2][1:])
    if t not in count:
        count[t] = 0
    
    with open(file_path, 'r') as file:
        for line in file:
            num = int(line.strip())
            if num > 10:
                count[t] += 1

i = 0
while i in count:
    print(count[i])
    i += 1
