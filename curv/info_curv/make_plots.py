import matplotlib.pyplot as plt
import numpy as np

def readfile(filename):
    f = open(filename,"r")
    len = []
    true_curv = []
    curv1 = []
    curv2 = []
    i = 0
    for line in f:
        if i==0:
            i+=1
            continue
        words = line.split()
        len.append(float(words[0]))
        true_curv.append(float(words[1]))
        curv1.append(float(words[2]))
        curv2.append(float(words[3]))
    return np.array(len), np.array(true_curv),np.array(curv1), np.array(curv2)

def make_plot_curv(len,true_curv, curv1, curv2, title,filename):
    plt.figure(figsize=(16,9))
    plt.plot(len,true_curv,linestyle = 'None',marker = 'x',label = "Истинная кривизна" )
    plt.plot(len,curv1, linestyle = "None",marker = 'x',label = "Обычные грани")
    plt.plot(len,curv2, linestyle = "None",marker = 'x', label = "Дуальные грани")
    plt.grid()
    plt.legend(fontsize = 16)
    plt.xlabel("len", fontsize = 16)
    plt.ylabel("Кривизна", fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.savefig(filename)

def make_plot_conv(nums_points, norms1,norms2, title,filename):
    plt.figure(figsize=(16,9))
    plt.plot(nums_points,norms1,linestyle = 'None',marker = 'x',label = "Отличие от истинной кривизны для обычных граней" )
    plt.plot(nums_points,norms2, linestyle = "None", marker = 'x',label = "Отличие от истинной кривизны для дуальных граней")
    plt.grid()
    plt.legend(fontsize = 16)
    plt.xlabel("количество точек", fontsize = 16)
    plt.ylabel("Норма", fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.savefig(filename)

def check_data(lens,true_curv,curv1,curv2):
    newlen = []
    newtrue = []
    newcurv1 = []
    newcurv2 = []
    for i in range(len(lens)):
        if(np.isnan(true_curv[i]) or np.isnan(curv1[i]) or np.isnan(curv2[i])):
            continue
        elif (true_curv[i]>2 or curv1[i]>2 or curv2[i]>2):
            continue
        elif (true_curv[i]<-2 or curv1[i]<-2 or curv2[i]<-2):
            continue
        else:
            newlen.append(lens[i])
            newtrue.append(true_curv[i])
            newcurv1.append(curv1[i])
            newcurv2.append(curv2[i])
    
    return newlen,newtrue,newcurv1, newcurv2

def all(fig,dist):
    filename1 = "./curv/" + fig + "s/" + fig + "_" + dist+ "_n"
    filename2 = "./pictures/" + fig + "s/" + fig + "_" + dist + "_n"
    isDelta = fig == 'delta'
    filename3 = "./pictures/" + fig + "s/" + fig + "_" + dist + "_conv.png"
    title = "Зависимость кривизны от длины"
    end1 = ".txt"
    end2 = ".png"
    indmax1 = 0
    indmax2 = 0
    norm1 = []
    norm2 = []
    num_points = np.array([i for i in range(5,101,5)])
    new_points = []
    when_plot = [5,10,20,25,30,80,100]
    j= 0
    for i in num_points:
        lens, true_curv, curv1,curv2 = readfile(filename1+str(i)+end1)
        if(isDelta):
            lens,true_curv,curv1,curv2 = check_data(lens,true_curv,curv1,curv2)
        new_points.append(i)
        indmax1 = np.argmax(np.abs(curv1 - true_curv))
        indmax2 = np.argmax(np.abs(curv1 - true_curv))
        norm1.append(curv1[indmax1] - true_curv[indmax1])
        norm2.append(curv2[indmax2] - true_curv[indmax2])
        if(j<len(when_plot)):
            if(when_plot[j] == i):
                j+=1
                make_plot_curv(lens,true_curv,curv1,curv2,title,filename2+str(i)+end2)
    make_plot_conv(num_points, norm1,norm2,"Зависимость ошибки от количества точек",filename3)

all("circle","uni")
all("circle","unipert")
#all("circle","gauss")
plt.close("all")

all("ellipse","uni")
all("ellipse","unipert")
#all("ellipse","gauss")
plt.close("all")

all("spiral","uni")
all("spiral","unipert")
#all("spiral","gauss")
plt.close("all")

#all("delta","uni")
#all("delta","unipert")
#all("delta","gauss")
#plt.close("all")

all("loop","uni")
all("loop","unipert")
#all("loop","gauss")
plt.close("all")

all("parabola","uni")
all("parabola","unipert")
#all("parabola","gauss")
plt.close("all")
