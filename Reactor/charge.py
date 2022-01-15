import sys

with open("temp.out", "r") as f:
    lines = f.readlines()
    for x in lines:
        if x.startswith('     #   Z          covCN'):
            num = lines.index(x)
            break
    charge = 0
    order = 5
    index1 = 0
    index2 = 0
    neighbors = [0,0,0]
    for i in range(num+1, num+41):
        line = ' '.join(lines[i].strip("\n").split())
        data = line.split(" ")
        if data[2] == "C":
            if float(data[4]) > charge:
                charge = float(data[4])
                index1 = data[0]
    #print(index)
    for i in range(num+53, num+110):
        line = ' '.join(lines[i].strip("\n").split())
        data = line.split(" ")
        try:
            if (data[2] == "C"):
                if (float(data[3]) < order):
                    order = float(data[3])
                    index2 = data[0]
        except:
            continue
    if (index1 == index2):
        print(index1)
    else:
        print("0")
