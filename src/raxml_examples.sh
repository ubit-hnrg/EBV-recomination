raxmlHPC -m GTRCAT -s /home/ariel/Projects/Gutierrez/EBV-recomb/JVI-resubmition/splited2/splited_1113-1192.fas -n TEST -p 12345 -b 12345 -N 10



raxmlHPC-PTHREADS -m GTRGAMMA -s /home/ariel/Projects/Gutierrez/EBV-recomb/JVI-resubmition/splited2/splited_1113-1192.fas -n 1113_n50 -p 12345 -N 50 -T 4 -f a -x 12345



############# 
#raxmlHPC-PTHREADS -m GTRGAMMA -s /home/ariel/Projects/Gutierrez/EBV-recomb/JVI-resubmition/splited2/splited_1113-119#2.fas -n whoo -p 12345 -N 20 -T 4 -f a -x 12345


raxmlHPC-PTHREADS -m GTRGAMMA -s /home/ariel/Projects/Gutierrez/EBV-recomb/JVI-resubmition/splited2/splited_1193-1335.fas -n test-1193 -p 12345 -N 20 -T 4 -x 12345 -f a 
raxmlHPC-PTHREADS -m GTRGAMMA -s /home/ariel/Projects/Gutierrez/EBV-recomb/JVI-resubmition/splited2/splited_1113-1192.fas -n test-1113 -p 12345 -N 20 -T 4 -x 12345 -f a  ## asi imprimom el arbol consenso



# el output que contiene el soporte para los nodos es: 
RAxML_bootstrap.XXX  , con XXX el nombre de la corrida (-n option)

#este file es un bunch de arboles, ergo se lo puedo pasar a otra corrida para calcular las metricas de incongruencia. 


# comparo (calculo IC y TC scores) 
raxmlHPC-PTHREADS -L MRE -m GTRGAMMA -n IC-test -T 4 -f i -t RAxML_bestTree.test-1113 -z RAxML_bootstrap.test-1193	

# 



