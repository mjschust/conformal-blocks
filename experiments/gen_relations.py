

rel_str = ""

for i in range(1, 6):
    rel_str += "y_" + str(i) + "^2 - 1, "

for i in range(1, 6):
    for j in range(i+1, 6):
        rel_str += "y_" + str(i) + "*y_" + str(j) +" - 2, "

for i in range(1, 6):
    for j in range(i+1, 6):
        for k in range(1, 6):
            if k == i or k == j: continue

            rel_str += "d_" + str(i) + str(j) + "*y_" + str(k) + " - 1, "

for i in range(1, 6):
    for j in range(i+1, 6):
        rel_str += "d_" + str(i) + str(j) + "^2 + 1, "

for i in range(1, 6):
    for j in range(i+1, 6):
        for k in range(1, 6):
            if k == i or k == j: continue
            for l in range(k+1, 6):
                if l == i or l == j: continue

                rel_str += "d_" + str(i) + str(j) + "*d_" + str(k) + str(l) +" - 1, "

print(rel_str[:-2])