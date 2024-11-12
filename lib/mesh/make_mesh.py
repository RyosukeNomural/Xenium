# Meshing of rectangular domain
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# default value
MESH_SIZE = 10
params=sys.argv
x_left = float(params[1])
x_right = float(params[2])
y_bottom = float(params[3])
y_top = float(params[4])

x_len = x_right - x_left # length in x-direction
y_len = y_top - y_bottom # length in y-direction
x_num = int((x_len + MESH_SIZE) // MESH_SIZE)  # division number in x-direction(切り上げ)
y_num = int((y_len + MESH_SIZE) // MESH_SIZE) # division number in y-direction
x0 = x_left # x-coordinate at left bottom of the domain
y0 = y_bottom # y-coordinate at left bottom of the domain

npoin = (x_num+1)*(y_num+1)
mesh_num = x_num*y_num
node = np.zeros([4,mesh_num],dtype=int)
x = np.zeros([2,npoin],dtype=float)


mesh_numbers = []
x_pointers, y_pointers = [], []
mesh_lefts, mesh_rights, mesh_bottoms, mesh_tops = [], [], [], [] # メッシュのxの小さい方,大きい方,y座標の小さい方,大きいほうをしまう
mesh_number = 0 # 何番目のメッシュに該当するか((0,0)→(0,1)→(0,2)の順番で1,2,3と振っていく)

for i in range(y_num - 1):
    bottom = y_bottom + i * MESH_SIZE #  # y座標の小さい方を10ずつ検索
    top = bottom + MESH_SIZE # y座標の大きい方を10ずつ検索
    y_pointer = int(bottom // MESH_SIZE) # y座標の小さい方を10で割った値をポインターにする

    for j in range(x_num - 1):
        left = x_left + j * MESH_SIZE
        right = left + MESH_SIZE
        x_pointer = int(left // MESH_SIZE) # x座標の小さい方を10で割った値をポインターにする
        mesh_number +=1

        mesh_bottoms.append(bottom)
        mesh_tops.append(top)
        mesh_lefts.append(left)
        mesh_rights.append(right)
        mesh_numbers.append(mesh_number)
        x_pointers.append(x_pointer)
        y_pointers.append(y_pointer)
        print(f"{mesh_number}:{left}_{right}_{bottom}_{top}")

df = pd.DataFrame(
        data={'mesh_pointer': mesh_numbers, 'left': mesh_lefts, 'right': mesh_rights, 
            'bottom':mesh_bottoms, 'top': mesh_tops, 'x_pointer': x_pointers, 'y_pointer': y_pointers},
        columns=['mesh_pointer', 'left', 'right', 'bottom', 'top', 'x_pointer', 'y_pointer'],
        index=None
    )

print(df)
df.to_csv('sample_mesh.csv', index=False)

