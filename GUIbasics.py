import sys
#from Tkinter import Tk, Label, Button, LEFT, RIGHT,W, StringVar, PhotoImage, Y
from Tkinter import *
class FirstGUI:
	LABEL_TEXT = [
		     "Default Values are already present",
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
		
		#self.greet_button = Button(master, text="Greet", command=self.greet)
		#self.greet_button.pack(side=LEFT)
		#self.greet_button.grid(row=5,column=0)
		
		#self.close_button = Button(master, text="Close", command=master.quit)
		#self.close_button.pack(side=RIGHT)
		#self.close_button.grid(row=5,column=1)
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


f = open(sys.argv[1],"w")

def msgentry():
	f.write("PHYSICAL_PROBLEM= "+prob.get()+"\n")
	f.write("KIND_TURB_MODEL= "+turb.get()+"\n")
	f.write("MATH_PROBLEM= DIRECT \n")
	f.write("RESTART_SOL= "+rest.get()+"\n")
	f.write("REGIME_TYPE= "+reg.get()+"\n")
	f.write("SYSTEM_MEASUREMENTS= "+si.get()+"\n")
	f.write("REYNOLDS_NUMBER= "+re.get()+"\n")
	f.write("MACH_NUMBER= "+ma.get()+"\n")
	f.write("AOA= "+aoa.get()+"\n")
	f.write("SIDESLIP_ANGLE= "+ss.get()+"\n")
	f.write("DISCARD_INFILES= "+infile.get()+"\n")
	f.write("FIXED_CL_MODE= "+clmode.get()+"\n")
	f.write("TARGET_CL= "+tcl.get()+"\n")
	f.write("INIT_OPTION= "+init.get()+"\n")
	f.write("FREESTREAM_OPTION= "+fs.get()+"\n")
	f.write("FREESTREAM_PRESSURE= "+prfs.get()+"\n")
	f.write("FREESTREAM_TEMPERATURE= "+tempfs.get()+"\n")
	f.write("REYNOLDS_LENGTH= "+relen.get()+"\n")
	f.write("FREESTREAM_DENSITY= "+denfs.get()+"\n")
	f.write("FREESTREAM_VELOCITY= ("+velx.get()+","+vely.get()+","+velz.get()+")\n")
	

PROB_OPTIONS = [
"EULER",
"NAVIER_STOKES",
"WAVE_EQUATION",
"HEAT_EQUATION",
"FEM_ELASTICITY",
"POISSON_EQUATION"
]

TURB_OPTIONS = [
"NONE",
"SA",
"SA_NEG",
"SST"
]

YESNO_OPTIONS = [
"NO",
"YES"
]

REGIME_OPTIONS = [
"COMPRESSIBLE",
"INCOMPRESSIBLE"
]

SI_OPTIONS = [
"SI",
"US"
]

INIT_OPTIONS = [
"REYNOLDS",
"TD_CONDITIONS"
]

FS_OPTIONS = [
"TEMPERATURE_FS",
"DENSITY_FS"
]

index = 1
prob = StringVar()
prob.set(PROB_OPTIONS[0])

turb = StringVar()
turb.set(TURB_OPTIONS[0])

rest = StringVar()
rest.set(YESNO_OPTIONS[0])

reg = StringVar()
reg.set(REGIME_OPTIONS[0])

si = StringVar()
si.set(SI_OPTIONS[0])

infile = StringVar()
infile.set(YESNO_OPTIONS[0])

clmode = StringVar()
clmode.set(YESNO_OPTIONS[0])

init = StringVar()
init.set(INIT_OPTIONS[0])

fs = StringVar()
fs.set(FS_OPTIONS[0])

re = StringVar()
re.set(3E6)
ma = StringVar()
ma.set(0.8)
aoa = StringVar()
aoa.set(1.25)
ss = StringVar()
ss.set(0.0)
tcl = StringVar()
tcl.set(0.8)
prfs = StringVar()
prfs.set(101325.0)
tempfs = StringVar()
tempfs.set(288.15)
relen = StringVar()
relen.set(1.0)
denfs = StringVar()
denfs.set(1.2886)
velx = StringVar()
velx.set(1.0)
vely = StringVar()
vely.set(0.0)
velz = StringVar()
velz.set(0.0)


label1 = Label(root, text="Problem Type: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,prob,*PROB_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="Turbulence Model: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,turb,*TURB_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="Restart Solution: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,rest,*YESNO_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="Regime Type: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,reg,*REGIME_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="System of Measurements: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,si,*SI_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label3 = Label(root, text="-------Compressible Freestream definition-------")
label3.grid(row=index,columnspan=3)
index = index+2

label2 = Label(root, text="Reynolds number: ")
label2.grid(row=index,column=0)
entry1 = Entry(root,textvariable=re)
entry1.grid(row=index,column=1)
index=index+1

label2 = Label(root, text="Mach number: ")
label2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=ma)
entry2.grid(row=index,column=1)
index=index+1

label2 = Label(root, text="Angle of Attack: ")
label2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=aoa)
entry2.grid(row=index,column=1)
index=index+1

label2 = Label(root, text="Sideslip angle: ")
label2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=ss)
entry2.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="Discard info files: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,infile,*YESNO_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="Fixed Cl mode: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,clmode,*YESNO_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label2 = Label(root, text="Target Cl: ")
label2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=tcl)
entry2.grid(row=index,column=1)
index=index+1

label3 = Label(root, text="Init option to choose between Reynolds and thermodynamic quantities")
label3.grid(row=index,columnspan=3)
index = index+1
label1 = Label(root, text="Intitialization Type: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,init,*INIT_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label3 = Label(root, text="Freestream option to choose between density and temperature for initialising solution")
label3.grid(row=index,columnspan=3)
index = index+1
label1 = Label(root, text="Freestream Option: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,fs,*FS_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label2 = Label(root, text="Freestream Pressure: ")
label2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=prfs)
entry2.grid(row=index,column=1)
index=index+1

label2 = Label(root, text="Freestream Temperature: ")
label2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=tempfs)
entry2.grid(row=index,column=1)
index=index+1

label2 = Label(root, text="Reynolds length: ")
label2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=relen)
entry2.grid(row=index,column=1)
index=index+1

label3 = Label(root, text="-------Incompressible Freestream definition-------")
label3.grid(row=index,columnspan=3)
index = index+2

label2 = Label(root, text="Freestream Density: ")
label2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=denfs)
entry2.grid(row=index,column=1)
index=index+1

label3 = Label(root, text="Velocity(x,y,z)")
label3.grid(row=index,column=0)
index = index+1

entry2 = Entry(root,textvariable=velx)
entry2.grid(row=index,column=0)
entry2 = Entry(root,textvariable=vely)
entry2.grid(row=index,column=1)
entry2 = Entry(root,textvariable=velz)
entry2.grid(row=index,column=2)
index=index+1

button = Button(root,text="Enter",command=msgentry)
button.grid(row=index,column=0)

button = Button(root,text="Exit",command=root.quit)
button.grid(row=index,column=1)


root.mainloop()
