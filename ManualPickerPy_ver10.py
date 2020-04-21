# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 18:17:19 2016

@author: tanakayumiko
@tuned by keisuke tsukada for python3 & PyQt5
"""

"""
**ウェルの数を自由に指定できるようにしたい
**グラフのプロットをもっと早くする
**時刻を取得しグラフに反映する
**回収したサンプルの情報を記録する
  どれをバックグラウンドにするかを選択する手間を省く
**解析したウェルの顕微鏡での番号が表示されるようにする（xmlファイルを読み込む必要）
  輝度を解析して分泌しているウェルを見つけるアルゴリズムを開発したい
**python3系に合わせて一部書き直し
"""
import sys
import csv
import numpy as np
from datetime import datetime
from PyQt5 import QtGui, QtCore ,QtWidgets
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import readFile_v3
import Analyze
import readXML
import datetime
import matplotlib.dates as mdates
#import sendGmail

class Run(QtWidgets.QWidget):
    global Roi_Data_file
    
    def __init__(self):
        super(Run, self).__init__()
        self.initUI()
        self.analyze()
        self.read()
            
    def initUI(self):
        self.analyze_btn = QtWidgets.QPushButton('Analyze scan:%d'%scan, self)
        self.analyze_btn.clicked.connect(self.update)
        self.read_btn = QtWidgets.QPushButton('Read', self)
        self.read_btn.clicked.connect(self.read)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.analyze_btn)
        self.layout.addWidget(self.read_btn)
        self.setLayout(self.layout)
        self.timer = QtCore.QTimer()
        self.timer.setInterval(100)#60*1000)
        #self.timer.timeout.connect(self.update)
        
    def analyze(self):
        global rawData
        global rawTime
        global pic
        global scan
        rawData, rawTime, pic, scan = Analyze.Analyze(scan, protein, nd2_file, Roi_Data_file, status_file, time_file, n)
        self.analyze_btn.setText('Analyze scan:%d'%scan)
        scan += 1

    def read(self):
        global Data
        global Status
        global Time
        global back
        global pic
        global timelag
        global Grid
        global scan
        global body
        global suggestion2
        try:
            Status.to_csv(status_file)
        except:
            pass
        Data, Status, Time = \
        readFile_v3.readFile(Roi_Data_file, status_file, time_file, rawData, rawTime)
        suggestion,Status = readFile_v3.back_fit_suggest(scan, Data, max(scan-15,1),Status)#readFile_v2.fitting(scan, Data, timelag, threthold)
        try:
            for i in range(n):
                for j in range(4):
                    Grid[i].Button[j].setStyleSheet('QPushButton {background-color: %s }'%color[int(Status.ix[i+1,j+1]['Status'])].name())
            """
            p =0
            fig = plt.figure(figsize=(4,3))
            for i in suggestion:
                if Status.ix[i[0],i[1]]['Status'] == 0.0 and i in suggestion2:
                    Grid[i[0]-1].Button[i[1]-1].setStyleSheet('QPushButton {background-color: #ff0000 }')
                    body = body + '%d_%d,\n'%(i[0],i[1])
                    p +=1
                    ax = fig.add_subplot(2, 2, p)
                    ax.set_title('%d_%d'%(i[0],i[1]))
                    ax.plot(Data.ix[i[0],i[1]][scan-20:scan], 'r-o', markersize=2)
                    #ax.plot(range(scan-3,scan),Data.ix[i[0],i[1]][scan-3:scan],'o')
                    ax.set_ylim([scan-20, scan])                    
                    ax.set_ylim([min(Data.ix[i[0],i[1]][scan-20:scan])-1, max(2,max(Data.ix[i[0],i[1]][scan-20:scan]))])
                    #plt.savefig('/Volumes/SINGLECELL/sampledata/20161007_max_test/plot.png')
            plt.close('all')
            """
            #if body != '分泌が疑われるサンプル：\n':
                #msg = sendGmail.create_message(ADDRESS, to_addr, subject, body, mime, attach_file) 
                #送信
                #sendGmail.send(ADDRESS, to_addr, msg)
                #body = '分泌が疑われるサンプル：\n'
        except:
            print('Error in read()')
        suggestion2 = suggestion       

    def update(self):
        try:
            self.analyze()
            self.read()
        except:
            print('Error in update()')
            
class MatplotlibWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(MatplotlibWidget, self).__init__(parent)
        
        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axis = self.figure.add_subplot(111)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.canvas)
        self.label = QtWidgets.QLabel(self)
        self.label.setGeometry(160, 40, 500, 30)
        self.label.setFont(QtGui.QFont('SanSerif', 30))
        
class Plot(QtWidgets.QWidget):
    def __init__(self):
        super(Plot, self).__init__()
        self.initUI()
        
    def initUI(self):
        self.plot_widget = MatplotlibWidget(self)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.plot_widget)

    @QtCore.pyqtSlot(int, int)
    def singlePlot(self, address, well):
        global Time
        global Status
        Time.columns = Time.columns.astype(int)
        if mode == 0:
            t = Time.loc[address]/60000
            t=t.sort_index()
            plot_len = min(len(t),20)
            plot_t = t[-plot_len:]
            plot_t = plot_t.values
            print(plot_t)
            line = Data.ix[address,well].values
            plot_line = line[-plot_len:]
            print(Data.ix[address,well])
            
            # Formatterでx軸の日付ラベルを月・日に設定
            #xfmt = mdates.DateFormatter("%H/%M")
            
            self.plot_widget.axis.clear()
            self.plot_widget.axis.grid()
            #self.plot_widget.axis.xaxis.set_major_formatter(xfmt)
            #self.plot_widget.axis.set_xlim([max(len(line)-20,0), len(line)+1])
            self.plot_widget.axis.set_ylim([min(min(plot_line),-0.5), max(max(plot_line),2)])
            slope,intercept,stderr = Status.loc[address,well]['slope'],\
            Status.loc[address,well]['intercept'],Status.loc[address,well]['stderr']
            self.plot_widget.axis.plot([plot_t[0], plot_t[-3]],\
            [slope*(max(scan-15,1))+intercept,slope*max(scan-3,0)+intercept],'--')
            self.plot_widget.axis.plot([plot_t[0], plot_t[-1]],\
            [slope*(max(scan-15,1))+intercept+stderr*sigma,slope*(scan)+intercept+stderr*sigma],'--')
            self.plot_widget.axis.plot([plot_t[-1]-30, plot_t[-1]-30],[min(min(plot_line),-0.5), max(max(plot_line),2)])
            self.plot_widget.axis.plot([plot_t[-1]-60, plot_t[-1]-60],[min(min(plot_line),-0.5), max(max(plot_line),2)])
            self.plot_widget.axis.plot([plot_t[-1]-90, plot_t[-1]-90],[min(min(plot_line),-0.5), max(max(plot_line),2)])

            self.plot_widget.axis.plot(plot_t, plot_line, '-o')#Data.ix[address,well],'-o')
            try:
                self.plot_widget.axis.plot([plot_t[-3],plot_t[-2],plot_t[-1]],\
                [plot_line[-3],plot_line[-2],plot_line[-1]],'o')#[(8/7)*(line[-3]/2+line[-4]/4+line[-5]/8),(8/7)*(line[-2]/2+line[-3]/4+line[-4]/8),(8/7)*(line[-1]/2+line[-2]/4+line[-3]/8)],'o')
            except:
                pass
            self.plot_widget.canvas.draw()
            self.plot_widget.label.setText('%d_%d : # %s'%(address, well, position[address-1][-3:]))
 
    @QtCore.pyqtSlot(int, int)
    def picture(self, address, well):
        pass
        
        if mode == 0:
            if   well == 1: a,b,c,d = 0,512,0,512
            elif well == 2: a,b,c,d = 0,512,510,1022
            elif well == 3: a,b,c,d = 512,1024,0,512
            elif well == 4: a,b,c,d = 512,1024,510,1022
            
            pic.default_coords['c'] = protein
            image = np.array(pic[address-1][a:b,c:d])
            self.plot_widget.axis.clear()
            self.plot_widget.axis.imshow(image)
            self.plot_widget.canvas.draw()
            self.plot_widget.label.setText('%d_%d : #%s'%(address, well, position[address-1][-3:]))
        
    @QtCore.pyqtSlot(int, int)
    def merge(self, address, well):
        pass
        
        if mode == 0:
            if   well == 1: a,b,c,d = 0,512,0,512
            elif well == 2: a,b,c,d = 0,512,510,1022
            elif well == 3: a,b,c,d = 512,1024,0,512
            elif well == 4: a,b,c,d = 512,1024,510,1022
            
            pic.default_coords['c'] = 0
            image = np.array(pic[address-1][a:b,c:d])
            self.plot_widget.axis.clear()
            self.plot_widget.axis.imshow(image, cmap = 'gray')
            self.plot_widget.canvas.draw()
            self.plot_widget.label.setText('%d_%d : #%s'%(address, well, position[address-1][-3:]))
        
#分泌検出アルゴリズムのパラメータ設定＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾＾
class Parameter(QtWidgets.QWidget):
    def __init__(self):
        super(Parameter, self).__init__()
        self.initUI()
        
    def initUI(self):
        #パラメーター3つ
        self.Difference = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)#輝度の差分
        self.Timelag = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)#いくつ前のタイムポイントと差分をとるか
        self.Threthold = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)#最初の撮影回との差分の閾値
        
        self.layout = QtWidgets.QVBoxLayout()
        
        #バーの最大値・最小値・初期値
        self.Difference.setMinimum(0)
        self.Difference.setMaximum(40)
        self.Difference.setValue(difference*10)
        
        self.Timelag.setMinimum(0)
        self.Timelag.setMaximum(10)
        self.Timelag.setValue(timelag)
        
        self.Threthold.setMinimum(0)
        self.Threthold.setMaximum(1000)
        self.Threthold.setValue(threthold*100)
        #ラベル
        self.Dlabel = QtWidgets.QLabel(self)
        self.Dlabel.setText('Difference %.1f'%difference)
        self.Tlabel = QtWidgets.QLabel(self)
        self.Tlabel.setText('Timelag %d'%timelag)
        self.Plabel = QtWidgets.QLabel(self)
        self.Plabel.setText('Threthold %.2f'%threthold)
        #バーとラベルの埋め込み
        self.layout.addWidget(self.Dlabel)
        self.layout.addWidget(self.Difference)
        self.layout.addWidget(self.Tlabel)
        self.layout.addWidget(self.Timelag)
        self.layout.addWidget(self.Plabel)
        self.layout.addWidget(self.Threthold)
        self.setLayout(self.layout)
        #バーの値を変更した時changeValueにとぶ
        self.Difference.valueChanged[int].connect(self.changeValue)
        self.Timelag.valueChanged[int].connect(self.changeValue)
        self.Threthold.valueChanged[int].connect(self.changeValue)
        
        self.show()
        
    def changeValue(self, value):
        global difference
        global timelag
        global threthold
        
        sender = self.sender()
         
        if sender == self.Difference:
            difference = value/10.0
            self.Dlabel.setText('Difference %.1f'%difference)
        elif sender == self.Timelag:
            timelag = value
            self.Tlabel.setText('Timelag %d'%timelag)
        elif sender == self.Threthold:
            threthold = value/100.0
            self.Plabel.setText('Threthold %.2f'%threthold)
            

class Address(QtWidgets.QWidget):
    trigger = QtCore.pyqtSignal(int, int)

    def __init__(self):
        super(Address, self).__init__()
        self.initUI()
        
    def initUI(self):
        global address
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        #左上（開始位置）の座標
        self.xstart = 1
        self.ystart = 1
        #ボタンのサイズ
        btnSize = 22
        #同じ番地内での仕切りのあつみ
        interval = 2
        #インスタンスが作られた時のaddressの値を自身のaddressに格納する
        self.address = address
        #各ボタンをクリックすると、ボタンの色が今選択中の色に変わる
        self.btn1 = QtWidgets.QPushButton(str(address) + '_1', self)
        self.btn2 = QtWidgets.QPushButton(str(address) + '_2', self)
        self.btn3 = QtWidgets.QPushButton(str(address) + '_3', self)
        self.btn4 = QtWidgets.QPushButton(str(address) + '_4', self)
        
        self.Button = [self.btn1, self.btn2, self.btn3, self.btn4]
        
        self.loc = [[self.xstart, self.ystart, btnSize, btnSize],\
                    [self.xstart + interval + btnSize, self.ystart, btnSize, btnSize],\
                    [self.xstart, self.ystart + interval + btnSize, btnSize, btnSize],\
                    [self.xstart + interval + btnSize, self.ystart + interval + btnSize, btnSize, btnSize]]
        for i in range(4):
            self.Button[i].setGeometry(self.loc[i][0], self.loc[i][1], self.loc[i][2], self.loc[i][3])
            self.Button[i].well = i+1
            self.Button[i].address = self.address
            self.Button[i].setStyleSheet\
            ('QPushButton {background-color: %s }'%color[int(Status.ix[self.address,i+1]['Status'])].name())
            self.Button[i].clicked.connect(self.singleplot)
            self.Button[i].clicked.connect(self.setStatus)
            
    def singleplot(self):
        sender = self.sender()
        self.trigger.emit(sender.address, sender.well)
        
    def setStatus(self):
        global color
        global status
        sender = self.sender()
        if mode == 1:
            (A, W) = (sender.address, sender.well)
            sender.setStyleSheet('QPushButton {background-color: %s }'%color[status].name())
            Status.ix[A,W] = status
            if status == 3:
                clock = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
                with open(sample_file, 'a') as f:
                    s_writer = csv.writer(f, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    s_writer.writerow([A]+[W]+[position[A-1]]+[clock])

class setStatus(QtWidgets.QWidget):
    def __init__(self):
        super(setStatus,self).__init__()
        self.initUI()
        
    def initUI(self):
        vbox = QtWidgets.QVBoxLayout(self)
        #setStatusModeにするボタン
        self.setStatusbutton = QtWidgets.QPushButton('setStatusMode', self)
        self.setStatusbutton.setCheckable(True)
        self.setStatusbutton.clicked[bool].connect(self.setMode)
        vbox.addWidget(self.setStatusbutton)
        #1細胞である
        self.singlecell = QtWidgets.QRadioButton('SingleCell', self)
        self.singlecell.toggled[bool].connect(self.selectColor)
        #self.singlecell.toggle()#はじめはここにチェックが入っている
        #QtCore.QObject.connect(self.singlecell, QtCore.SIGNAL('toggled(bool)'), self.selectColor)
        vbox.addWidget(self.singlecell)

        #ごみ入りである
        self.omit = QtWidgets.QRadioButton('Omit', self)
        self.omit.toggled[bool].connect(self.selectColor)
        #QtCore.QObject.connect(self.omit, QtCore.SIGNAL('toggled(bool)'), self.selectColor)
        vbox.addWidget(self.omit)

        #バックグラウンドにすべきである
        self.background = QtWidgets.QRadioButton('BackGround', self)
        self.background.toggled[bool].connect(self.selectColor)
        #QtCore.QObject.connect(self.background, QtCore.SIGNAL('toggled(bool)'), self.selectColor)
        vbox.addWidget(self.background)
        
        #回収済み
        self.collected = QtWidgets.QRadioButton('Collected', self)
        self.collected.toggled[bool].connect(self.selectColor)
        #QtCore.QObject.connect(self.collected, QtCore.SIGNAL('toggled(bool)'), self.selectColor)
        vbox.addWidget(self.collected)
        
        self.setLayout(vbox)
        self.show()
        
    def setMode(self, pressed):
        global mode
        if pressed:
            mode = 1
        else:
            mode = 0
    #ラジオボタンを押すたびに選択中の色と文字列が変わる 
    def selectColor(self):
        global status
        sender = self.sender()
        if sender.text() == 'SingleCell':
            status = 0
        elif sender.text() == 'Omit':
            status = 1
        elif sender.text() == 'BackGround':
            status = 2
        elif sender.text() == 'Collected':
            status = 3
            
def main():
    #str,画像を解析し平均輝度を求めた結果が書いてあるcsvファイルのパス
    global Roi_Data_file
    #str,顕微鏡から吐き出されるnd2ファイルのパスの接頭語（このうしろにつづく番号はあとで付け足す）
    global nd2_file
    #str,XY position（ND aquisition)撮影画面の番号、通し番号じゃない#がついてる方
    global position
    #str,ピックした時刻が記録されるcsvファイルのパス
    global time_file
    #str,ピックしたウェルの場所が記録されるcsvファイルのパス
    global sample_file
    #str,各ウェルが、「1細胞、2:omit、3:空、4:回収済み」のどの状態なのかが記録されるcsvファイルのパス
    global status_file
    #int,XY position (ND aquisition)通し番号の方
    global address
    #DataFrame,各ウェルの輝度平均の時系列データ
    global Data
    #list,ボタンの色
    global color
    #int(0/1),mode=1のときウェルボタンをつつけばウェルのstatus変更、mode=0のときならウェルの輝度プロット表示
    global mode
    #int,各撮影回で得られる総画面数。実験回により異なる
    global n
    #int,現在までの撮影回の数
    global scan
    #int,ウェルが「1細胞、2:omit、3:空、4:回収済み」のどの状態にいるかを示す
    global status
    #float,今の平均輝度とtimelag回前の撮影回の平均輝度との差がこれ以上だと分泌を始めたとしてアラートする
    global difference
    #int,
    global timelag
    #float,今の平均輝度と最初の撮影回の平均輝度との差がこれ以上だと分泌を始めたとしてアラートする
    global threthold
    #list,n個のAddressオブジェクトの羅列
    global Grid
    #int,何番目のフィルターの画像を解析するか(DIA,CF660R,Cy5のとき、Cy5を解析するなら2)
    global protein
    #str,メールの宛先
    global to_addr
    #str,メールタイトル
    global subject
    #str,メール本文
    global body
    #int,直近10回前後の撮影回の輝度を1次関数にフィッティングし、そのstd*sigmaを3回超えるとアラートする（たしか）
    global sigma
    #str,送信元アドレス（gmailのセキュリティーの設定を下げないと送れない)
    global ADDRESS
    #dic,メールに画像を添付するときの設定
    global mime
    #dic,メールに添付するファイルの情報（名前、パス）
    global attach_file
    
    app = QtWidgets.QApplication(sys.argv)
    """
    Roi_Data_file = r'E:\sampledata_ROI_Data.csv'
    position_file = r'Z:/RealTimeExport/20161007_Quad2_hILC2_pick3.xml'
    nd2_file = r'Z:/images/Image'
    time_file = r'C:/Users/YTanaka/Desktop/20161007/time.csv'
    sample_file = r'C:/Users/YTanaka/Desktop/20161007/sample.csv'
    status_file = r'C:/Users/YTanaka/Desktop/20161007/status.csv'
    """
    filepath = r'C:/Users/gomak/OneDrive - The University of Tokyo/Documents/sampledata_'
    Roi_Data_file = filepath + r'\ROI_Data.csv'
    position_file = filepath + r'\20161007_Quad2_hILC2_pick3.xml'
    nd2_file = filepath + r'\nd2file\Image'#'
    time_file = filepath + r'\time.csv'
    sample_file = filepath + r'\sample.csv'
    status_file = filepath + r'\status.csv'
    
    n = 159#変更する
    protein = 2#フィルター
    scan = 27#新しくGUIを立ち上げる時は変更
    pos_n = 996
    #件名と本文
    subject = '2号機より'
    body = '分泌が疑われるサンプル：\n'
    #宛先アドレス
    to_addr = "gomakanabun@gmail.com" 
    #送信元アドレス
    ADDRESS = 'komamehanamuguri@gmail.com'
    #添付ファイル設定(text.txtファイルを添付)
    mime={'type':'image', 'subtype':'png'}
    attach_file={'name':'plot.png', 'path':'/Volumes/SINGLECELL/sampledata/20161007_max_test/plot.png'}
 
    #メッセージの作成(添付ファイルあり)
    mode = 0
    status = 2
    difference = 1
    timelag = 0
    threthold = 1
    sigma = 50
    color = [QtGui.QColor(100, 30, 30),#赤
            QtGui.QColor(100, 100, 100),#灰色
            QtGui.QColor(255, 255, 200),#クリーム
            QtGui.QColor(200, 255, 255)]#水色
            
    position = readXML.readXML(position_file, pos_n)
    
    main_panel = QtWidgets.QWidget()
    main_panel.resize(10000,1000)
    
    panel = QtWidgets.QWidget()
    panel2 = QtWidgets.QWidget()
    read_widget = Run()
    status_widget = setStatus()
    graph_widget = Plot()
    para_widget = Parameter()
    pic_widget = Plot()
    merge_widget = Plot()
    tab_widget = QtWidgets.QTabWidget()
    well_widget = QtWidgets.QWidget()
    
    main_panel_layout = QtWidgets.QHBoxLayout()
    panel_layout = QtWidgets.QVBoxLayout()
    panel2_layout = QtWidgets.QHBoxLayout()
    tab_layout = QtWidgets.QHBoxLayout()
    well_layout = QtWidgets.QGridLayout()

    #番地番号をまず付け足していき、そのまま各番地のインスタンスを作成
    Grid = []

    for i in np.arange(1, n+1):
        address = i
        #まずは番地番号をAddressの最後尾につけたす
        Grid.append(str(i))
        #Addressの最後尾にWellのインスタンスを作成
        Grid[-1] = Address()
        Grid[-1].trigger.connect(graph_widget.singlePlot)
        Grid[-1].trigger.connect(pic_widget.picture)
        Grid[-1].trigger.connect(merge_widget.merge)

        #その番地番号に従い、グリット状にウィジェットを配置（左上から右下へ）
        well_layout.addWidget(Grid[-1], (i-1)//10, (i-1)%10)

    tab_widget.addTab(graph_widget, 'graph')
    tab_widget.addTab(pic_widget, 'picture')
    tab_widget.addTab(merge_widget, 'merge')
    tab_layout.addWidget(tab_widget)
    panel_layout.addWidget(read_widget)
    panel2_layout.addWidget(status_widget)
    panel2_layout.addWidget(para_widget)
    panel2.setLayout(panel2_layout)
    panel_layout.addWidget(panel2)
    panel_layout.addWidget(tab_widget)
    panel.setLayout(panel_layout)
    well_widget.setLayout(well_layout)
    main_panel_layout.addWidget(well_widget)
    main_panel_layout.addWidget(panel)
    main_panel.setLayout(main_panel_layout)
    main_panel.show()
    read_widget.timer.start()
    sys.exit(app.exec_())
    
if __name__=='__main__':
    main()