from numpy import *
from tqdm import *
from random import randint
from matplotlib import *
from tkinter import *
import numpy as np
import matplotlib.pyplot as plt
import pyautogui
import os

#to do : bouton exit dessin

class First_layers: # fct d'activation : Relu
    def __init__(self, weights, biases):
        self.weights = weights
        self.biases = biases

    def forward(self,inputs):   
        self.inputs = inputs
        self.output_not_activated = np.dot(self.weights, inputs) + self.biases
        self.output = np.maximum(0, self.output_not_activated)

    def backward(self, input_gradient, learning_rate):
        input_gradient = input_gradient * deriv_relu(self.output_not_activated)
        self.output_gradient = np.dot(self.weights.T, input_gradient)
        self.weights -= np.dot(input_gradient, self.inputs.T) * learning_rate * 1/np.shape(input_gradient)[1]
        self.biases -= np.mean(input_gradient, axis = 1, keepdims = True) * learning_rate 

class Last_layer: # fct d'activation : Softmax
    def __init__(self, weights, biases): 
        self.weights = weights
        self.biases = biases

    def forward(self,inputs): 
        self.inputs = inputs
        output_not_activated = np.dot(self.weights, inputs) + self.biases
        exp_values = np.exp(output_not_activated - np.max(output_not_activated, axis = 0, keepdims = True))
        self.output = exp_values / np.sum(exp_values, axis = 0, keepdims = True) 
    
    def backward(self, y_pred, y_true, learning_rate): # fct de perte : Categorical-crossentropy
        input_gradient = y_pred - y_true
        self.output_gradient = np.dot(self.weights.T, input_gradient)
        self.weights -= np.dot(input_gradient, self.inputs.T) * learning_rate * 1/np.shape(input_gradient)[1]
        self.biases -= np.mean(input_gradient, axis = 1, keepdims = True) * learning_rate

def open():
    global nn, acc_input, acc_label, dstest
    weights = np.load("weights.npy", [], True)
    biases = np.load("biases.npy", [], True)
    dstest = np.load("test.npy", [], True)
    n = len(weights)
    acc_input = np.array([dstest[i][1] for i in range(10000)]).T
    acc_label = np.array([dstest[i][0] for i in range(10000)]).T
    nn = [First_layers(weights[j], biases[j]) for j in range(n-1)] + [Last_layer(weights[n-1], biases[n-1])]

def create(npc): # npc : neurones par couche
    n = len(npc)
    info = np.array([npc, 0, [0]], dtype = object)
    biases = np.array([np.random.randn(npc[j+1],1) for j in range(n-1)], dtype = object)
    weights = np.array([0.1 * np.random.randn(npc[j+1],npc[j]) for j in range(n-1)], dtype = object)
    np.save("info.npy", info)
    np.save("weights.npy", weights)
    np.save("biases.npy", biases)

def get_info():
    info = np.load("info.npy", [], True)
    print(info[0])
    plt.plot(range(len(info[2])), info[2]) 
    plt.title("best accuracy : "+str(info[1]))
    plt.show()

def save():
    try:
        info = np.load("info.npy", [], True) 
        folder = str(round(info[1]*100,2))+" "+str(info[0])
        os.makedirs(folder)
        os.rename("biases.npy",folder+"/biases.npy")
        os.rename("weights.npy",folder+"/weights.npy")
        os.rename("info.npy",folder+"/info.npy")
    except OSError:
        raise FileExistsError("pas de fichiers à sauvergarder")

def delete():
    try:
        os.remove("biases.npy")
        os.remove("weights.npy")
        os.remove("info.npy")
    except OSError:
        raise FileExistsError("pas de fichiers à supprimer")

def load(n):
    folders = os.listdir()
    folders.sort()
    folder = folders[n]
    os.rename(folder+"/biases.npy", "biases.npy")
    os.rename(folder+"/weights.npy", "weights.npy")
    os.rename(folder+"/info.npy", "info.npy")

def resave(n, rename_folder = False):
    info = np.load("info.npy", [], True)
    folders = os.listdir()
    folders.sort()
    folder = folders[n]
    if len(os.listdir(folder)) == 0:
        if rename_folder:
            new_folder = str(round(info[1]*100,2))+" "+str(info[0])
            os.rename(folder, new_folder)
            folder = new_folder
        os.rename("biases.npy",folder+"/biases.npy")
        os.rename("weights.npy",folder+"/weights.npy")
        os.rename("info.npy",folder+"/info.npy")
    else:
        raise FileExistsError("des fichiers se trouvent déjà dans ce dossier")

def ls():
    for file in os.listdir():
        if file[0] in "12345667890":
            print(file)

def train(learning_rate, t_samples, nb, nbbt = 1):
    dstrain = np.load("train.npy", [], True)
    info = np.load("info.npy", [], True)
    best_acc = info[1]
    open()
    for i in tqdm(range(nb)):
        for j in range(nbbt):
            choix = [randint(0,59999) for i in range(t_samples)] 
            inputs = np.array([dstrain[choix[i]][1] for i in range(t_samples)]).T
            labels = np.array([[i==dstrain[choix[j]][0] for i in range(10)] for j in range(t_samples)]).T
            output = forpropagation(inputs)
            backpropagation(output, labels, learning_rate)
        acc = accuracy()
        if acc == 0 or best_acc - acc > 0.3:
            raise ValueError("l'accuracy a chuté, diminuer le learning rate pourrait régler le problème")
        info[2].append(acc)
        if acc>best_acc:
            best_acc = acc
            best_weights = np.array([(nn[j].weights).copy() for j in range(len(nn))], dtype = object)
            best_biases = np.array([(nn[j].biases).copy() for j in range(len(nn))], dtype = object)
    info[1] = best_acc
    np.save("info.npy", info)
    try :
        np.save("weights.npy", best_weights)
        np.save("biases.npy", best_biases)
    except UnboundLocalError:
        raise ValueError("pas de meilleure accuracy trouvé")
    
def test(nb, only_false = False): 
    open()
    for i in range(nb):
        choix = randint(0,9999)
        input = np.array([dstest[choix][1]]).T
        label = dstest[choix][0]
        output = forpropagation(input)
        if not only_false or np.argmax(output)!=label:
            image = np.reshape(input,(28,28))
            plt.gray()
            plt.imshow(image)
            plt.title("prediction : "+str(np.argmax(output))+" label : "+str(label))
            plt.show(block=False)
            plt.pause(10)
    plt.close()

def accuracy():
    output = forpropagation(acc_input)
    return np.mean(np.argmax(output, axis = 0) == acc_label)

def accuracy1():
    for i in tqdm(range(1)):
        output = forpropagation(acc_input)
        a = np.mean(np.argmax(output, axis = 0) == acc_label)
    return a

def accuracy2():
    for i in tqdm(range(1)):
        input = np.array([dstest[i][1] for i in range(10000)]).T
        label = np.array([dstest[i][0] for i in range(10000)]).T
        output = forpropagation(input)
        a = np.mean(np.argmax(output, axis = 0) == label)
    return a

def accuracy3():
    for i in tqdm(range(1)):
        l = []
        for i in range(10000):
            input = np.array([dstest[i][1]]).T
            label = np.array(dstest[i][0]).T
            output = forpropagation(input)
            l.append(np.max(output)==output[label][0])
        a = np.mean(l)
    return a

def backpropagation(y_pred, y_true, learning_rate): 
    nn[-1].backward(y_pred, y_true, learning_rate)
    for i in range(len(nn)-2,-1,-1):
        nn[i].backward(nn[i+1].output_gradient, learning_rate)
    
def forpropagation(input):
    nn[0].forward(input)
    for i in range(1,len(nn)):
        nn[i].forward(nn[i-1].output)
    return nn[-1].output

def deriv_relu(x):
    return x>0

def dessin():
    global mode, condition
    open()
    root = Tk()
    draw = [[0 for i in range(28)] for i in range(28)]
    input = [[0 for i in range(28)] for i in range(28)]
    canvas = Canvas(root, background="black", width=700, height=800)
    text = canvas.create_text((200,745), text = str([0]*10) + " none", fill = "white")
    for i in range(0,701,25):
        canvas.create_line((0, i), (700, i), fill='grey25', width=1)
        canvas.create_line((i, 0), (i, 700), fill='grey25', width=1)
    root.focus_force()
    canvas.pack()
    condition = False
    mode = "pencil"

    def pointeur():
        if condition:
            xy = [xycan[0] + (pyautogui.position()[0] - xyscr[0]), xycan[1] + (pyautogui.position()[1]-xyscr[1])]
            if (0 < xy[0] - xy[0]%25 + 25 < 701 and 0 < xy[1] - xy[1]%25 + 25 < 701) and (mode == "eraser" or (mode == "pencil" and input[xy[1]//25][xy[0]//25] != 255)):
                n = 2
                for i in range(-n,n+1):
                    for j in range(-n,n+1):
                        if sqrt(i**2 + j**2)<=n:
                            if 0 < xy[0] - xy[0]%25 + 25 + i*25 < 701 and 0 < xy[1] - xy[1]%25 + 25 + j*25 < 701:
                                if mode == "pencil":
                                    c = min(int((1 - (sqrt(i**2 + j**2)/n)) * 255) + input[xy[1]//25 + j][xy[0]//25 + i], 255)
                                    rgb = "#%02x%02x%02x" % (c,c,c)
                                    if c != 0:
                                        input[xy[1]//25 + j][xy[0]//25 + i] = c
                                        if draw[xy[1]//25 + j][xy[0]//25 + i] == 0:
                                            draw[xy[1]//25 + j][xy[0]//25 + i] = canvas.create_rectangle((xy[0] - xy[0]%25 + i*25, xy[1] - xy[1]%25 + j*25), (xy[0] - xy[0]%25+25 + i*25, xy[1] - xy[1]%25 + 25 + j*25), fill = rgb, outline = rgb)
                                        else:
                                            canvas.itemconfig(draw[xy[1]//25 + j][xy[0]//25 + i], fill = rgb)
                                            canvas.itemconfig(draw[xy[1]//25 + j][xy[0]//25 + i], outline = rgb)
                                if mode == "eraser" and draw[xy[1]//25 + j][xy[0]//25 + i] != 0:
                                    canvas.delete(draw[xy[1]//25 + j][xy[0]//25 + i])
                                    draw[xy[1]//25 + j][xy[0]//25 + i]=0
                                    input[xy[1]//25 + j][xy[0]//25 + i]=0          
                input1 = np.reshape(np.array(input), (784,1))
                pred = forpropagation(input1)
                canvas.itemconfig(text, text = str(np.ndarray.tolist(np.reshape(np.round(pred, 1),(10,))))+ " " + str(np.argmax(pred)))
        root.after(1, pointeur)

    def beg(event):
        global xycan,condition,xyscr
        xycan,condition,xyscr = [event.x,event.y],True,pyautogui.position()

    def end(event):
        global condition
        condition = False

    def pencil():
        global mode
        bg2.configure(bg = "white")
        bg1.configure(bg = "green")
        mode = "pencil"

    def eraser():
        global mode
        bg2.configure(bg = "green")
        bg1.configure(bg = "white")
        mode = "eraser"

    def reset():
        for i in range(28):
            for j in range(28):
                if draw[i][j] != 0:
                    canvas.delete(draw[i][j])
                    draw[i][j] = 0
                    input[i][j] = 0
        canvas.itemconfig(text, text = str([0]*10) + " none")

    bg1 = LabelFrame(root, bd=6, bg="green")
    bg1.place(x=400,y=28*25+15)
    icon1 = PhotoImage(file='pencil.png')
    Button1 = Button(bg1, height=50, width=50, image=icon1, command=pencil)
    Button1.pack()

    bg2 = LabelFrame(root, bd=6, bg="white")
    bg2.place(x=500, y=28*25+15)
    icon2 = PhotoImage(file='eraser.png')
    Button2 = Button(bg2, height=50, width=50, image=icon2, command=eraser)
    Button2.pack()

    bg3 = LabelFrame(root, bd=6, bg="white")
    bg3.place(x=600, y=28*25+22)
    Button3 = Button(bg3, height=1, width=5, text="delete", command=reset)
    Button3.pack()

    pointeur()
    canvas.bind("<Button-1>", beg)
    canvas.bind("<ButtonRelease-1>", end)
    root.mainloop()

