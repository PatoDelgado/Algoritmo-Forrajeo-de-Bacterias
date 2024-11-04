from bacteria import bacteria
from chemiotaxis import chemiotaxis

import numpy

poblacion = []
path = "/Users/patriciodelgado/Desktop/BFOA-main/multiFasta.fasta"
numeroDeBacterias = 10
numRandomBacteria = 1
iteraciones = 50
tumbo = 1                                              #numero de gaps a insertar 
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)                #mejor bacteria   
tempBacteria = bacteria(path)            #bacteria temporal para validaciones
original = bacteria(path)                #bacteria original sin gaps
globalNFE = 0      #numero de evaluaciones de la funcion objetivo

dAttr= 0.15 #0.1
wAttr= 0.25 #0.2
hRep=0.12
wRep= 12    #10


def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction
    
def validaSecuencias(path, veryBest):
    #clona a veryBest en tempBacteria   
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    #descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-","")
    #tempBacteria.tumboNado(1)    

    #valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return
            #nueva funcion 
def nuevo_movimiento_bacterias(bacteria, tumbo, nado):
    print(f"Fitness antes del movimiento: {bacteria.fitness}")
    bacteria.tumboNado(tumbo)  # movimiento tumbo como en el original
    bacteria.autoEvalua()  # evaluación del fitness
    print(f"Fitness después del movimiento: {bacteria.fitness}")

    if bacteria.fitness < 0.5:
        print("Fitness bajo, aplicando movimiento amplio (doble nado)")
        bacteria.tumboNado(nado * 2)  # si el fitness es bajo mover más
    else:
        print("Fitness alto, aplicando movimiento preciso (mitad tumbo)")
        bacteria.tumboNado(tumbo // 2)  # si el fitness es alto movimientos finos

    print(f"Secuencia de la bacteria después del movimiento: {bacteria.matrix.seqs}")


for i in range(numeroDeBacterias):                                            #poblacion inicial
    poblacion.append(bacteria(path))


for _ in range(iteraciones):                                                  #numero de iteraciones  
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)
        #bacteria.tumboNado(nado)
        bacteria.autoEvalua()  
        #print("blosumScore: ",bacteria.blosumScore)
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)                 #d_attr, w_attr, h_rep, w_rep):
    globalNFE += chemio.parcialNFE 
    best = max(poblacion, key=lambda x: x.fitness)
    if (veryBest == None) or (best.fitness > veryBest.fitness):
        clonaBest(veryBest, best)
    print("interaccion: ",veryBest.interaction,"fitness: ",veryBest.fitness, " NFE:",globalNFE )
    
    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)                #inserta  bacterias aleatorias
    print("poblacion: ",len(poblacion))


veryBest.showGenome()
validaSecuencias(path, veryBest)