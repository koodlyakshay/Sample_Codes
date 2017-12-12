import sys
#from Tkinter import Tk, Label, Button, LEFT, RIGHT,W, StringVar, PhotoImage, Y
from Tkinter import *
class FirstGUI:
	LABEL_TEXT = [
		     "This is our first GUI",
		     "Second line",
		     "Third line",
		     "Fourth line",
		     "Fifth line",
		]
	def __init__(self, master):
				
		self.master = master
		master.title("A simple GUI")
			
		#self.label = Label(master, text="Hello World")
		#self.label.pack()
		self.label_index = 0
		self.label_text = StringVar()
		self.label_text.set(self.LABEL_TEXT[self.label_index])
		self.label = Label(master, textvariable=self.label_text)
		self.label.bind("<Button-1>", self.cycle_label_text)
		self.label.pack()
		
		self.label.grid(columnspan=3, sticky=W)
		
		self.greet_button = Button(master, text="Greet", command=self.greet)
		#self.greet_button.pack(side=LEFT)
		self.greet_button.grid(row=4,column=0)
		
		self.close_button = Button(master, text="Close", command=master.quit)
		#self.close_button.pack(side=RIGHT)
		self.close_button.grid(row=4,column=1)
	def greet(self):
	    print("Greetings")
	 
	def cycle_label_text(self, event):
		self.label_index += 1
		self.label_index %=len(self.LABEL_TEXT[self.label_index])
		self.label_text.set(self.LABEL_TEXT[self.label_index])
def printmsg():
	print("HAha")
	
	
root = Tk()
my_gui = FirstGUI(root)		


menu_bar = Menu(root)

file_menu = Menu(menu_bar, tearoff=0)
file_menu.add_command(label="Quit", command=root.quit)
file_menu.add_command(label="Exit", command=root.quit)
file_menu.add_command(label="HAha", command=printmsg)

menu_bar.add_cascade(label="File", menu=file_menu)

root.config(menu=menu_bar)

re = StringVar()
pr = StringVar()
f = open("test2.txt","w")

def msgentry():
	f.write("REYNOLDS NUMBER = "+re.get()+"\n")
	f.write("PRANDTL NUMBER = "+pr.get()+"\n")



label1 = Label(root, text="Reynolds number: ")
label1.grid(row=0,column=0)
entry1 = Entry(root,textvariable=re)
entry1.grid(row=0,column=1)

label2 = Label(root, text="Prandtl number: ")
label2.grid(row=1,column=0)
entry2 = Entry(root,textvariable=pr)
entry2.grid(row=1,column=1)






button = Button(root,text="Enter",command=msgentry)
button.grid(row=2,column=1)


root.mainloop()
