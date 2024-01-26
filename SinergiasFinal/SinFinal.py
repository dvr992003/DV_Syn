import sys
import os
import numpy as np
import time
import matplotlib.pyplot as plt
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog
from PyQt5 import uic

class Application(QDialog):
    def __init__(self):

       
        self.data_matrix = [[]] 
        self.chan = []          
        self.V_matrix = []      
        self.H_matrix = [[]]    
        self.W_matrix = [[]]    
        self.NNMF_rank = 0      
        self.iter_limit = 1     
        self.error_thr = 0.001 
        self.error_mode = True  

        self.fig = []
        

        super(Application,self).__init__()
        uic.loadUi("gui.ui",self)

        self.chan_selector_tm = QtWidgets.QMenu(self)

        self.browse_btn.clicked.connect(self.browse_files)
        self.load_btn.clicked.connect(self.load_data)
        self.rank_selector.valueChanged.connect(self.rank_selected)
        self.chan_selector.view().pressed.connect(self.channel_selected)
        self.init_btn.clicked.connect(self.init_test)
        self.display_btn.clicked.connect(self.show_results)
        self.save_btn.clicked.connect(self.save_results)
        

    def browse_files(self):                                     

        fname = QFileDialog.getOpenFileName(self,'Open File', 'C:', 'TXT files (*.txt);;CSV Files ( *.csv)')
        self.filename_text.setText(fname[0])                    
        if(fname[0] != ''):
            self.load_btn.setEnabled(True)                     
        else:
            self.load_btn.setEnabled(False)                     

        self.rank_selector.setEnabled(False)                   
        self.iter_selector.setEnabled(False)                    
        self.thr_selector.setEnabled(False)                     
        self.thr_ticker.setEnabled(False)                      
        self.display_btn.setEnabled(False)                      
        self.save_btn.setEnabled(False)                        
        self.init_btn.setEnabled(False)                       

        self.size_selector.setEnabled(False)
        self.size_selector.setValue(0)

        self.status_label.setText("Idle")                       

    def load_data(self):                                        

       
        CSVData = open(self.filename_text.text())
        self.data_matrix = np.loadtxt(CSVData, delimiter=",")   
        self.data_indicator.setText("Data has " + str(len(self.data_matrix)) + " channels and " + str(len(self.data_matrix[0])) + " samples.")

        
        self.NNMF_rank = 0                                      
        self.rank_selector.setEnabled(True)                     
        self.rank_selector.setMaximum(0)                        

        self.size_selector.setEnabled(True)
        self.size_selector.setMaximum(len(self.data_matrix[0]))
        self.size_selector.setValue(len(self.data_matrix[0]))

        self.iter_selector.setEnabled(True)                     
        self.thr_selector.setEnabled(True)                     
        self.thr_ticker.setEnabled(True)                        
        self.display_btn.setEnabled(False)                      
        self.save_btn.setEnabled(False)                         

        self.status_label.setText("Idle")                      

        self.chan_selector_tm.clear()                          
        self.chan_selector.clear()                             
        self.chan_selector.addItem('--Select Channels--')      

        for i in range(len(self.data_matrix)):                   
            self.chan_selector.addItem('Channel %s' % i)
            item = self.chan_selector.model().item(i+1, 0)
            item.setCheckState(QtCore.Qt.Unchecked)
            action = self.chan_selector_tm.addAction('Channel %s' % i)
            action.setCheckable(True)
    
    def rank_selected(self):                                 
        self.NNMF_rank = self.rank_selector.value()
        self.init_btn.setEnabled(self.NNMF_rank > 0)           


    def channel_selected(self, index):                         
        item = self.chan_selector.model().itemFromIndex(index)  
        if(item.row() != 0):                                   
            if item.checkState() == QtCore.Qt.Checked:          
                item.setCheckState(QtCore.Qt.Unchecked)
                self.rank_selector.setMaximum(self.rank_selector.maximum()-1)
            else:
                item.setCheckState(QtCore.Qt.Checked)
                self.rank_selector.setMaximum(self.rank_selector.maximum()+1)
        
    def init_test(self):                                       
        self.chan = [int for i in range(len(self.data_matrix))] 
        self.V_matrix = []

        for i in range(len(self.data_matrix)):                 
            self.chan[i] = self.chan_selector.model().item(i+1, 0).checkState() == QtCore.Qt.Checked
        
        chan_list = []

        for i in range(len(self.data_matrix)):                 
            if(self.chan[i]):
                chan_list.append(i)
                self.V_matrix.append(list(self.data_matrix[i][0:(self.size_selector.value()-1)]))

        self.H_matrix = np.random.rand(self.NNMF_rank,len(self.V_matrix[0]))    
        self.W_matrix = np.random.rand(len(self.V_matrix),self.NNMF_rank)       
        self.iter_limit = self.iter_selector.value()                           
        self.error_thr = self.thr_selector.value()                              
        self.error_mode = self.thr_ticker.isChecked()                           
        self.display_btn.setEnabled(False)                                     
        self.save_btn.setEnabled(False)                                         
        self.H_matrix, self.W_matrix, self.fig = calculate(self.V_matrix, self.W_matrix, self.H_matrix, self.iter_limit,self.error_thr,self.error_mode,self.status_label,chan_list)

        self.display_btn.setEnabled(True)                                       
        self.save_btn.setEnabled(True)                                          

    def show_results(self): 
        for f in self.fig:
            f.show()

    def save_results(self): 

        filepath = os.path.dirname(self.filename_text.text())
        basename = os.path.basename(self.filename_text.text())

        print("a")
        print(filepath + "H_" + basename)

        np.savetxt(filepath + "\\H_" + basename, self.H_matrix, delimiter=',')
        np.savetxt(filepath + "\\W_" + basename, self.W_matrix, delimiter=',')
    
        
app = QApplication(sys.argv)
nnmf_app = Application()
widget=QtWidgets.QStackedWidget()
widget.addWidget(nnmf_app)
widget.setFixedWidth(400)
widget.setFixedHeight(300)
widget.show()
nnmf_app.show()

def mult_martix(m1, m2):

    c = [[0 for col in range(len(m2[0]))] for row in range(len(m1))]

    for i in range(len(m1)):
        for j in range(len(m2[0])):
            for k in range(len(m2)):
                c[i][j] += m1[i][k] * m2[k][j]
                
    return c

def div_martix(m1, m2):

    c = [[0 for col in range(len(m2[0]))] for row in range(len(m1))]

    for i in range(len(m1)):
        for j in range(len(m1[0])):
            c[i][j] = m1[i][j] / m2[i][j]
                
    return c

def tran_matrix(m1):

    c = [[0 for col in range(len(m1))] for row in range(len(m1[0]))]

    for i in range(len(m1)):
        for j in range(len(m1[0])):
            c[j][i] = m1[i][j]

    return c

def dot_product(m1, m2):

    c = [[0 for col in range(len(m2[0]))] for row in range(len(m1))]

    for i in range(len(m1)):
        for j in range(len(m1[0])):
            c[i][j] = m1[i][j] * m2[i][j]

    return c

def add_matrix(m1,m2,add_mode):

    c = [[0 for col in range(len(m2[0]))] for row in range(len(m1))]

    for i in range(len(m1)):
        for j in range(len(m1[0])):
            c[i][j] = m1[i][j] + m2[i][j]*(-1+2*add_mode)

    return c

def thr_matrix(m1,threshold,cummulative):

    valid = False

    if(cummulative):
        c = 0
        for i in range(len(m1)):
            for j in range(len(m1[0])):
                c += abs(m1[i][j])
        if(c <= threshold):
            valid = True
    
    else:
        valid = True
        c = 0
        for i in range(len(m1)):
            for j in range(len(m1[0])):
                c = abs(m1[i][j])
                if(c > threshold):
                    valid = False
    
    return valid


def calculate(V, W, H, max_iter, threshold, cummulative, label, chan_list):

    found_result = False
    num_iters = 0
    in_time = time.time()
    for n in range(max_iter):
        label.setText("Calculating iteration: " + str(n) + "/" + str(max_iter)) 
        app.processEvents()                                                     

        H_mult_num = mult_martix(tran_matrix(W),V)
        H_mult_den = mult_martix(mult_martix(tran_matrix(W),W),H)
        H_mult = div_martix(H_mult_num, H_mult_den)
        H = dot_product(H,H_mult)

        W_mult_num = mult_martix(V,tran_matrix(H))
        W_mult_den = mult_martix(mult_martix(W,H),tran_matrix(H))
        W_mult = div_martix(W_mult_num,W_mult_den)
        W = dot_product(W,W_mult)

        V_final = mult_martix(W,H)
        diff = add_matrix(V,V_final,False)

        if(thr_matrix(diff,threshold,cummulative)):
            found_result = True
            num_iters = n
            break

    out_time = time.time()
    if(found_result):  
        label.setText("Found in " + str(out_time-in_time) + "s after " + str(num_iters) + " passes.")
    else:
        label.setText("Not found in " + str(out_time-in_time) + "s after " + str(max_iter) + " passes.")

    f = [0]*len(V_final)

    for chan in range(len(V_final)):    
        f[chan] = plt.figure(chan)
        plt.plot(V[chan],label="L1")
        plt.plot(V_final[chan], label="L2")
        plt.plot(diff[chan], label="L3")
        plt.legend(('Original V','Reconstructed V','Error'))
        plt.title(f"Channel {chan_list[chan]}")
        plt.ylim(-1, 1)

    return H,W,f

sys.exit(app.exec_())

