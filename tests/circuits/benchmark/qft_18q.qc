version 1.0

qubits 18

# init qubits to |0...01> (set your own input state here)
.init
   x q[0]

# main circuit
.qft
   h q[0]
   cr q[1],q[0]
   cr q[2],q[0]
   cr q[3],q[0]
   cr q[4],q[0]
   cr q[5],q[0]
   cr q[6],q[0]
   cr q[7],q[0]
   cr q[8],q[0]
   cr q[9],q[0]
   cr q[10],q[0]
   cr q[11],q[0]
   cr q[12],q[0]
   cr q[13],q[0]
   cr q[14],q[0]
   cr q[15],q[0]
   cr q[16],q[0]
   cr q[17],q[0]
   h q[1]
   cr q[2],q[1]
   cr q[3],q[1]
   cr q[4],q[1]
   cr q[5],q[1]
   cr q[6],q[1]
   cr q[7],q[1]
   cr q[8],q[1]
   cr q[9],q[1]
   cr q[10],q[1]
   cr q[11],q[1]
   cr q[12],q[1]
   cr q[13],q[1]
   cr q[14],q[1]
   cr q[15],q[1]
   cr q[16],q[1]
   cr q[17],q[1]
   h q[2]
   cr q[3],q[2]
   cr q[4],q[2]
   cr q[5],q[2]
   cr q[6],q[2]
   cr q[7],q[2]
   cr q[8],q[2]
   cr q[9],q[2]
   cr q[10],q[2]
   cr q[11],q[2]
   cr q[12],q[2]
   cr q[13],q[2]
   cr q[14],q[2]
   cr q[15],q[2]
   cr q[16],q[2]
   cr q[17],q[2]
   h q[3]
   cr q[4],q[3]
   cr q[5],q[3]
   cr q[6],q[3]
   cr q[7],q[3]
   cr q[8],q[3]
   cr q[9],q[3]
   cr q[10],q[3]
   cr q[11],q[3]
   cr q[12],q[3]
   cr q[13],q[3]
   cr q[14],q[3]
   cr q[15],q[3]
   cr q[16],q[3]
   cr q[17],q[3]
   h q[4]
   cr q[5],q[4]
   cr q[6],q[4]
   cr q[7],q[4]
   cr q[8],q[4]
   cr q[9],q[4]
   cr q[10],q[4]
   cr q[11],q[4]
   cr q[12],q[4]
   cr q[13],q[4]
   cr q[14],q[4]
   cr q[15],q[4]
   cr q[16],q[4]
   cr q[17],q[4]
   h q[5]
   cr q[6],q[5]
   cr q[7],q[5]
   cr q[8],q[5]
   cr q[9],q[5]
   cr q[10],q[5]
   cr q[11],q[5]
   cr q[12],q[5]
   cr q[13],q[5]
   cr q[14],q[5]
   cr q[15],q[5]
   cr q[16],q[5]
   cr q[17],q[5]
   h q[6]
   cr q[7],q[6]
   cr q[8],q[6]
   cr q[9],q[6]
   cr q[10],q[6]
   cr q[11],q[6]
   cr q[12],q[6]
   cr q[13],q[6]
   cr q[14],q[6]
   cr q[15],q[6]
   cr q[16],q[6]
   cr q[17],q[6]
   h q[7]
   cr q[8],q[7]
   cr q[9],q[7]
   cr q[10],q[7]
   cr q[11],q[7]
   cr q[12],q[7]
   cr q[13],q[7]
   cr q[14],q[7]
   cr q[15],q[7]
   cr q[16],q[7]
   cr q[17],q[7]
   h q[8]
   cr q[9],q[8]
   cr q[10],q[8]
   cr q[11],q[8]
   cr q[12],q[8]
   cr q[13],q[8]
   cr q[14],q[8]
   cr q[15],q[8]
   cr q[16],q[8]
   cr q[17],q[8]
   h q[9]
   cr q[10],q[9]
   cr q[11],q[9]
   cr q[12],q[9]
   cr q[13],q[9]
   cr q[14],q[9]
   cr q[15],q[9]
   cr q[16],q[9]
   cr q[17],q[9]
   h q[10]
   cr q[11],q[10]
   cr q[12],q[10]
   cr q[13],q[10]
   cr q[14],q[10]
   cr q[15],q[10]
   cr q[16],q[10]
   cr q[17],q[10]
   h q[11]
   cr q[12],q[11]
   cr q[13],q[11]
   cr q[14],q[11]
   cr q[15],q[11]
   cr q[16],q[11]
   cr q[17],q[11]
   h q[12]
   cr q[13],q[12]
   cr q[14],q[12]
   cr q[15],q[12]
   cr q[16],q[12]
   cr q[17],q[12]
   h q[13]
   cr q[14],q[13]
   cr q[15],q[13]
   cr q[16],q[13]
   cr q[17],q[13]
   h q[14]
   cr q[15],q[14]
   cr q[16],q[14]
   cr q[17],q[14]
   h q[15]
   cr q[16],q[15]
   cr q[17],q[15]
   h q[16]
   cr q[17],q[16]
   h q[17]

# reverse qubit order
#.reverse
   swap q[0],q[17]
   swap q[1],q[16]
   swap q[2],q[15]
   swap q[3],q[14]
   swap q[4],q[13]
   swap q[5],q[12]
   swap q[6],q[11]
   swap q[7],q[10]
   swap q[8],q[9]

# statistics: 
#  total number of gates: 198
#    |- number of H    gates: 18
#    |- number of Rk   gates: 153
#    |- number of swap gates: 9(27 cnot gates)

