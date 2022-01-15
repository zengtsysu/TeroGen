import sys
input = sys.argv[1]
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
            #print(str(data[0]) + "  " + str(data[3]))
            '''
            if (data[0] == index) and (float(data[3]) < 3.9):
                neighbors = [data[5], data[8], data[11]]
                neighbors.sort()
                with open("cation.dat", "a") as f:
                    f.write("state"+str(input) + " " + str(index) + " " + str(neighbors[0]) + " " + str(neighbors[1]) + " " + str(neighbors[2]) + "\n")
                break
            '''
            if (data[2] == "C"):
                if (float(data[3]) < order):
                    order = float(data[3])
                    index2 = data[0]
        except:
            continue
    if (index1 == index2):
        with open("cation.dat", "a") as f:
            f.write("state"+str(input) + " " + str(index1) + "\n")
