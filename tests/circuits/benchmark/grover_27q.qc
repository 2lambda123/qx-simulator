version 1.0

qubits 27

.init
   x q[13]
   h q[0]
   h q[1]
   h q[2]
   h q[3]
   h q[4]
   h q[5]
   h q[6]
   h q[7]
   h q[8]
   h q[9]
   h q[10]
   h q[11]
   h q[12]
   h q[13]

# iterating 71 times
.grover(71)

   # oracle_128

   x q[0]
   x q[1]
   x q[2]
   x q[3]
   x q[4]
   x q[5]
   x q[6]
   x q[8]
   x q[9]
   x q[10]
   x q[11]
   x q[12]
   toffoli q[0],q[1],q[14]
   toffoli q[1],q[14],q[15]
   toffoli q[2],q[15],q[16]
   toffoli q[3],q[16],q[17]
   toffoli q[4],q[17],q[18]
   toffoli q[5],q[18],q[19]
   toffoli q[6],q[19],q[20]
   toffoli q[7],q[20],q[21]
   toffoli q[8],q[21],q[22]
   toffoli q[9],q[22],q[23]
   toffoli q[10],q[23],q[24]
   toffoli q[11],q[24],q[25]
   toffoli q[12],q[25],q[26]
   cnot q[26],q[13]
   toffoli q[12],q[25],q[26]
   toffoli q[11],q[24],q[25]
   toffoli q[10],q[23],q[24]
   toffoli q[9],q[22],q[23]
   toffoli q[8],q[21],q[22]
   toffoli q[7],q[20],q[21]
   toffoli q[6],q[19],q[20]
   toffoli q[5],q[18],q[19]
   toffoli q[4],q[17],q[18]
   toffoli q[3],q[16],q[17]
   toffoli q[2],q[15],q[16]
   toffoli q[1],q[14],q[15]
   toffoli q[0],q[1],q[14]
   x q[0]
   x q[1]
   x q[2]
   x q[3]
   x q[4]
   x q[5]
   x q[6]
   x q[8]
   x q[9]
   x q[10]
   x q[11]
   x q[12]
   h q[0]
   h q[1]
   h q[2]
   h q[3]
   h q[4]
   h q[5]
   h q[6]
   h q[7]
   h q[8]
   h q[9]
   h q[10]
   h q[11]
   h q[12]

   # inversion

   x q[0]
   x q[1]
   x q[2]
   x q[3]
   x q[4]
   x q[5]
   x q[6]
   x q[7]
   x q[8]
   x q[9]
   x q[10]
   x q[11]
   x q[12]
   h q[12]
   toffoli q[0],q[1],q[14]
   toffoli q[1],q[14],q[15]
   toffoli q[2],q[15],q[16]
   toffoli q[3],q[16],q[17]
   toffoli q[4],q[17],q[18]
   toffoli q[5],q[18],q[19]
   toffoli q[6],q[19],q[20]
   toffoli q[7],q[20],q[21]
   toffoli q[8],q[21],q[22]
   toffoli q[9],q[22],q[23]
   toffoli q[10],q[23],q[24]
   toffoli q[11],q[24],q[25]
   cnot q[25],q[12]
   toffoli q[11],q[24],q[25]
   toffoli q[10],q[23],q[24]
   toffoli q[9],q[22],q[23]
   toffoli q[8],q[21],q[22]
   toffoli q[7],q[20],q[21]
   toffoli q[6],q[19],q[20]
   toffoli q[5],q[18],q[19]
   toffoli q[4],q[17],q[18]
   toffoli q[3],q[16],q[17]
   toffoli q[2],q[15],q[16]
   toffoli q[1],q[14],q[15]
   toffoli q[0],q[1],q[14]
   h q[12]
   x q[0]
   x q[1]
   x q[2]
   x q[3]
   x q[4]
   x q[5]
   x q[6]
   x q[7]
   x q[8]
   x q[9]
   x q[10]
   x q[11]
   x q[12]
   h q[0]
   h q[1]
   h q[2]
   h q[3]
   h q[4]
   h q[5]
   h q[6]
   h q[7]
   h q[8]
   h q[9]
   h q[10]
   h q[11]
   h q[12]

.measure
   h q[13]
   measure q[13]
   display
