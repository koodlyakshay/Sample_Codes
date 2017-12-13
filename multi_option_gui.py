import sys
from Tkinter import *

def msgentry_comp():
	f.write("REGIME_TYPE= COMPRESSIBLE\n")
	f.write("REYNOLDS_NUMBER= "+re.get()+"\n")
	f.write("MACH_NUMBER= "+ma.get()+"\n")
	f.write("AOA= "+aoa.get()+"\n")
	f.write("SIDESLIP_ANGLE= "+ss.get()+"\n")
	f.write("DISCARD_INFILES= "+infile.get()+"\n")
	f.write("INIT_OPTION= "+init.get()+"\n")
	f.write("FREESTREAM_OPTION= "+fs.get()+"\n")
	f.write("FREESTREAM_PRESSURE= "+prfs.get()+"\n")
	f.write("FREESTREAM_TEMPERATURE= "+tempfs.get()+"\n")
	f.write("REYNOLDS_LENGTH= "+relen.get()+"\n")
	
def msgentry_stdair():
	f.write("FLUID_MODEL= STANDARD_AIR\n")
	f.write("GAMMA= 1.4\n")
	f.write("CRITICAL_TEMPERATURE= 131.0\n")
	f.write("GAS_CONSTANT= 287.058\n")
	f.write("CRITICAL_PRESSURE= 3588550.0\n")
	f.write("ACCENTRIC_FACTOR= 0.035\n")


def msgentry_flmodel():
	f.write("FLUID_MODEL= "+flopt.get()+"\n")
	f.write("GAMMA= "+gamma.get()+"\n")
	f.write("CRITICAL_TEMPERATURE= "+t_c.get()+"\n")
	f.write("GAS_CONSTANT= "+gas_cons.get()+"\n")
	f.write("CRITICAL_PRESSURE= "+p_c.get()+"\n")
	f.write("ACCENTRIC_FACTOR= "+a_f.get()+"\n")
	
def msgentry_viscmodel():
	f.write("VISCOSITY_MODEL= SUTHERLAND\n")
	f.write("MU_REF= "+mu_ref.get()+"\n")
	f.write("MU_T_REF= "+mu_tref.get()+"\n")
	f.write("SUTHERLAND_CONSTANT= "+sut_cons.get()+"\n")
	
def msgentry_consvisc():
	f.write("VISCOSITY_MODEL= CONSTANT_VISCOSITY\n")
	f.write("MU_CONS= "+mu_cons.get()+"\n")



def msgentry_incomp():
	f.write("REGIME_TYPE= INCOMPRESSIBLE\n")
	f.write("FREESTREAM_DENSITY= "+denfs.get()+"\n")
	f.write("FREESTREAM_VELOCITY= ("+velx.get()+","+vely.get()+","+velz.get()+")\n")
	f.write("FREESTREAM_VISCOSITY= "+visc.get()+"\n")

def msgentry_basic():
	f.write("PHYSICAL_PROBLEM= "+prob.get()+"\n")
	f.write("MATH_PROBLEM= DIRECT \n")
	f.write("RESTART_SOL= "+rest.get()+"\n")
	f.write("SYSTEM_MEASUREMENTS= "+si.get()+"\n")
	
def msgentry_refval():
	f.write("REF_NONDIMENSIONALISATION= "+refnondim.get()+"\n")
	f.write("REF_LENGTH= "+reflen.get()+"\n")
	f.write("REF_AREA= "+refar.get()+"\n")
	f.write("REF_ORIGIN_MOMENT_X="+refmomx.get()+"\n")
	f.write("REF_ORIGIN_MOMENT_Y="+refmomy.get()+"\n")
	f.write("REF_ORIGIN_MOMENT_Z="+refmomz.get()+"\n")
	f.write("FREESTREAM_VISCOSITY= "+visc.get()+"\n")
	
def msgentry_unsteady():
	f.write("UNSTEADY_SIMULATION= "+unst.get()+"\n")
	f.write("UNST_TIMESTEP= "+unst_delt.get()+"\n")
	f.write("UNST_TIME= "+unst_time.get()+"\n")
	f.write("UNST_CFL_NUMBER= "+unst_cfl.get()+"\n")
	f.write("UNST_INT_ITER= "+n_intiter.get()+"\n")
	f.write("UNST_RESTART_ITER= "+unst_rest.get()+"\n")
	
def msgentry_marker():
	f.write("MARKER_EULER= ("+euler.get()+")\n")
	f.write("MARKER_SYM= ("+sym.get()+")\n")
	f.write("MARKER_FAR= ("+far.get()+")\n")
	f.write("MARKER_INTERNAL= ("+internal.get()+")\n")
	f.write("MARKER_NEARFIELD= ("+nearfield.get()+")\n")
	f.write("MARKER_PRESSURE= ("+pressure.get()+")\n")
	f.write("MARKER_DIRICHLET= ("+dirichlet.get()+")\n")
	f.write("MARKER_NEUMANN= ("+neumann.get()+")\n")
	f.write("MARKER_OUTLET= ("+outlet.get()+","+backpr_val.get()+")\n")
	f.write("MARKER_HEATFLUX= ("+heatflux.get()+","+hflux_val.get()+")\n")
	f.write("MARKER_ISOTHERMAL= ("+isotherm.get()+","+isot_val.get()+")\n")
	
def msgentry_iden():
	f.write("MARKER_PLOTTING= ("+plotting.get()+")\n")
	f.write("MARKER_MONITORING= ("+monitor.get()+")\n")
	f.write("MARKER_DESIGNING= ("+evaluate.get()+")\n")
	

root = Tk()

f=open(sys.argv[1],"w")

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

FLUID_OPTIONS = [
"STANDARD_AIR",
"IDEAL_GAS",
"VW_GAS",
"PR_GAS"
]

NONDIM_OPTIONS = [
"DIMENSIONAL", 
"FREESTREAM_PRESS_EQ_ONE",
"FREESTREAM_VEL_EQ_MACH," 
"FREESTREAM_VEL_EQ_ONE"
]

UNST_OPTIONS = [
"NO", 
"TIME_STEPPING",
"DUAL_TIME_STEPPING-1ST_ORDER", 
"DUAL_TIME_STEPPING-2ND_ORDER",
"HARMONIC_BALANCE"
]

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

flopt = StringVar()
flopt.set(FLUID_OPTIONS[0])

refnondim = StringVar()
refnondim.set(NONDIM_OPTIONS[0])

unst = StringVar()
unst.set(UNST_OPTIONS[0])

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
visc = StringVar()
visc.set(1.853E-5)
dcldal = StringVar()
dcldal.set(0.0)
alphaupd = StringVar()
alphaupd.set(5)
reflen = StringVar()
reflen.set(1.0)
refar = StringVar()
refar.set(0.0)
sspan = StringVar()
sspan.set(1.0)
refmomx = StringVar()
refmomx.set(0.25)
refmomy = StringVar()
refmomy.set(0.0)
refmomz = StringVar()
refmomz.set(0.0)
gamma = StringVar()
gamma.set(1.4)
gas_cons = StringVar()
gas_cons.set(287.0)
t_c = StringVar()
t_c.set(131.0)
p_c = StringVar()
p_c.set(3588550.0)
a_f = StringVar()
a_f.set(0.035)
mu_cons = StringVar()
mu_cons.set(1.716E-5)
mu_ref = StringVar()
mu_ref.set(1.716E-5)
mu_tref = StringVar()
mu_tref.set(273.15)
sut_cons = StringVar()
sut_cons.set(110.4)
unst_delt = StringVar()
unst_delt.set(0.0)
unst_time = StringVar()
unst_time.set(50.0)
unst_cfl = StringVar()
unst_cfl.set(0.0)
n_intiter = StringVar()
n_intiter.set(200)
unst_rest = StringVar()
unst_rest.set(0)
euler = StringVar()
euler.set(NONE)
far = StringVar()
far.set(NONE)
sym = StringVar()
sym.set(NONE)
pressure = StringVar()
pressure.set(NONE)
neumann = StringVar()
neumann.set(NONE)
dirichlet = StringVar()
dirichlet.set(NONE)
internal = StringVar()
internal.set(NONE)
nearfield = StringVar()
nearfield.set(NONE)
interface = StringVar()
interface.set(NONE)
isotherm = StringVar()
isotherm.set(NONE)
isot_val = StringVar()
isot_val.set(0.0)
heatflux = StringVar()
heatflux.set(NONE)
hflux_val = StringVar()
hflux_val.set(0.0)
outlet = StringVar()
outlet.set(NONE)
backpr_val = StringVar()
backpr_val.set(0.0)
plotting = StringVar()
plotting.set(NONE)
monitor = StringVar()
monitor.set(NONE)
evaluate = StringVar()
evaluate.set(NONE)

def incompressible_window():
	incomp_win=Toplevel()
	incomp_win.title("Incompressible Freestream definition")
	index = 1

	label2 = Label(incomp_win, text="Freestream Density: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(incomp_win,textvariable=denfs)
	entry2.grid(row=index,column=1)
	index=index+1

	label3 = Label(incomp_win, text="Velocity(x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(incomp_win,textvariable=velx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(incomp_win,textvariable=vely)
	entry2.grid(row=index,column=1)
	entry2 = Entry(incomp_win,textvariable=velz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label2 = Label(incomp_win, text="Freestream viscosity: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(incomp_win,textvariable=visc)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(incomp_win,text="Enter values",command=msgentry_incomp)
	button.grid(row=index,column=0)
	button = Button(incomp_win,text="Back",command=incomp_win.destroy)
	button.grid(row=index,column=1)

def compressible_window():
	second_win = Toplevel()
	second_win.title('Compressible Options')
	index =1 
	label2 = Label(second_win, text="Reynolds number: ")
	label2.grid(row=index,column=0)
	entry1 = Entry(second_win,textvariable=re)
	entry1.grid(row=index,column=1)
	index=index+1

	label2 = Label(second_win, text="Mach number: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(second_win,textvariable=ma)
	entry2.grid(row=index,column=1)
	index=index+1

	label2 = Label(second_win, text="Angle of Attack: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(second_win,textvariable=aoa)
	entry2.grid(row=index,column=1)
	index=index+1

	label2 = Label(second_win, text="Sideslip angle: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(second_win,textvariable=ss)
	entry2.grid(row=index,column=1)
	index=index+1

	label1 = Label(second_win, text="Discard info files: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(second_win,infile,*YESNO_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1

	label3 = Label(second_win, text="Init option to choose between Reynolds and thermodynamic quantities")
	label3.grid(row=index,columnspan=3)
	index = index+1
	label1 = Label(second_win, text="Intitialization Type: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(second_win,init,*INIT_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1

	label3 = Label(second_win, text="Freestream option to choose between density and temperature for initialising solution")
	label3.grid(row=index,columnspan=3)
	index = index+1
	label1 = Label(second_win, text="Freestream Option: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(second_win,fs,*FS_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1

	label2 = Label(second_win, text="Freestream Pressure: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(second_win,textvariable=prfs)
	entry2.grid(row=index,column=1)
	index=index+1

	label2 = Label(second_win, text="Freestream Temperature: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(second_win,textvariable=tempfs)
	entry2.grid(row=index,column=1)
	index=index+1

	label2 = Label(second_win, text="Reynolds length: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(second_win,textvariable=relen)
	entry2.grid(row=index,column=1)
	index=index+1
    
    
	label1 = Label(second_win, text="Fluid Model: ")
	label1.grid(row=index,column=0)
	button = Button(second_win,text="STANDARD_AIR",command=msgentry_stdair)
	button.grid(row=index,column=1)
	button = Button(second_win,text="OTHER",command=flmodel_window)
	button.grid(row=index,column=2)
	index=index+1 
	
	label1 = Label(second_win, text="Viscosity Model: ")
	label1.grid(row=index,column=0)
	button = Button(second_win,text="SUTHERLAND",command=viscmodel_window)
	button.grid(row=index,column=1)
	button = Button(second_win,text="CONSTANT_VISCOSITY",command=get_consvisc)
	button.grid(row=index,column=2)
	index=index+1 
     
	button = Button(second_win,text="Enter values",command=msgentry_comp)
	button.grid(row=index,column=0)
	button = Button(second_win,text="Back",command=second_win.destroy)
	button.grid(row=index,column=1)
	
def get_consvisc():
	consvisc_win = Toplevel()
	consvisc_win.title("Constant Viscosity")
	index = 1
	label2 = Label(consvisc_win, text="Constant Viscosity: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(consvisc_win,textvariable=mu_cons)
	entry2.grid(row=index,column=1)
	index=index+1
	button = Button(consvisc_win,text="Write",command=msgentry_consvisc)
	button.grid(row=index,column=0)
	button = Button(consvisc_win,text="Back",command=consvisc_win.destroy)
	button.grid(row=index,column=1)
	
def viscmodel_window():
	viscmodel_win = Toplevel()
	viscmodel_win.title("Viscous Model definition")
	
	index = 1
	
	label2 = Label(viscmodel_win, text="Sutherland Viscosity Reference: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(viscmodel_win,textvariable=mu_ref)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(viscmodel_win, text="Sutherland Temperature Ref: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(viscmodel_win,textvariable=mu_tref)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(viscmodel_win, text="Sutherland constant: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(viscmodel_win,textvariable=sut_cons)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(viscmodel_win,text="Enter values",command=msgentry_viscmodel)
	button.grid(row=index,column=0)
	button = Button(viscmodel_win,text="Back",command=viscmodel_win.destroy)
	button.grid(row=index,column=1)



def flmodel_window():
	flmodel_win = Toplevel()
	flmodel_win.title("FLUID MODEL DEFINITION")
	
	index = 1
	
	label1 = Label(flmodel_win, text="Fluid Model: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(flmodel_win,flopt,*FLUID_OPTIONS)
	opt1.grid(row=index,column=1)
	index+=1
	
	label2 = Label(flmodel_win, text="Gamma: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(flmodel_win,textvariable=gamma)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(flmodel_win, text="Gas Constant: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(flmodel_win,textvariable=gas_cons)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(flmodel_win, text="Critical Temperature: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(flmodel_win,textvariable=t_c)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(flmodel_win, text="Critical Pressure: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(flmodel_win,textvariable=p_c)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(flmodel_win, text="Accentric Factor ")
	label2.grid(row=index,column=0)
	entry2 = Entry(flmodel_win,textvariable=a_f)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(flmodel_win,text="Enter values",command=msgentry_flmodel)
	button.grid(row=index,column=0)
	button = Button(flmodel_win,text="Back",command=flmodel_win.destroy)
	button.grid(row=index,column=1)
	
def define_refvalues():
	ref_win = Toplevel()
	ref_win.title("Marker definition")
	
	index=1
	
	label2 = Label(ref_win, text="Reference length: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(ref_win,textvariable=reflen)
	entry2.grid(row=index,column=1)
	index=index+1

	label3 = Label(ref_win, text="Reference Origin for moment(x,y,z)f")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(ref_win,textvariable=refmomx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(ref_win,textvariable=refmomy)
	entry2.grid(row=index,column=1)
	entry2 = Entry(ref_win,textvariable=refmomz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label2 = Label(ref_win, text="Reference Area: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(ref_win,textvariable=refar)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(ref_win, text="Aircraft Semi Span: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(ref_win,textvariable=sspan)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label1 = Label(ref_win, text="Non dimensionalisation Type: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(ref_win,refnondim,*NONDIM_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	button = Button(ref_win,text="Enter values",command=msgentry_refval)
	button.grid(row=index,column=0)
	button = Button(ref_win,text="Back",command=ref_win.destroy)
	button.grid(row=index,column=1)
	
def define_unsteady():
	unsteady_win = Toplevel()
	unsteady_win.title("UNSTEADY PARAMETERS DEFINITION")
	
	index = 1
	
	label1 = Label(unsteady_win, text="Type of unsteady Simulation: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(unsteady_win,unst,*UNST_OPTIONS)
	opt1.grid(row=index,column=1)
	index+=1
	
	label2 = Label(unsteady_win, text="Time step for dual time stepping: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(unsteady_win,textvariable=unst_delt)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(unsteady_win, text="Total physical time: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(unsteady_win,textvariable=unst_time)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(unsteady_win, text="Unsteady CFL for finest grid: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(unsteady_win,textvariable=unst_cfl)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(unsteady_win, text="Internal iterations: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(unsteady_win,textvariable=n_intiter)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(unsteady_win, text="Start restart at iteration: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(unsteady_win,textvariable=unst_rest)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(unsteady_win,text="Enter values",command=msgentry_unsteady)
	button.grid(row=index,column=0)
	button = Button(unsteady_win,text="Back",command=unsteady_win.destroy)
	button.grid(row=index,column=1)
	
def define_markers():
	marker_win = Toplevel()
	marker_win.title("Marker definition")
	
	index = 1
		
	label2 = Label(marker_win, text="Euler Wall: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=euler)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Farfield: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=far)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Symmetry: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=sym)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Internal: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=internal)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Nearfield: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=nearfield)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Interface: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=interface)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Pressure marker: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=pressure)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Neumann: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=neumann)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Dirichlet: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=dirichlet)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Heatflux(Name, Flux): ")
	label2.grid(row=index,column=0)
	index+=1
	entry2 = Entry(marker_win,textvariable=heatflux)
	entry2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=hflux_val)
	entry2.grid(row=index,column=1)
	index+=1
	
	label2 = Label(marker_win, text="Isothermal(Name, Temp(K)): ")
	label2.grid(row=index,column=0)
	index+=1
	entry2 = Entry(marker_win,textvariable=isotherm)
	entry2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=isot_val)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(marker_win, text="Outlet(Name, Back pressure(static)): ")
	label2.grid(row=index,column=0)
	index+=1
	entry2 = Entry(marker_win,textvariable=outlet)
	entry2.grid(row=index,column=0)
	entry2 = Entry(marker_win,textvariable=backpr_val)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(marker_win,text="Enter values",command=msgentry_marker)
	button.grid(row=index,column=0)
	button = Button(marker_win,text="Back",command=marker_win.destroy)
	button.grid(row=index,column=1)
	
def define_iden_markers():
	iden_win = Toplevel()
	iden_win.title("Marker definition")
	
	index = 1
		
	label2 = Label(iden_win, text="Markers in suface flow file: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(iden_win,textvariable=plotting)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(iden_win, text="Markers to evalate non-dim coefficients: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(iden_win,textvariable=monitor)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(iden_win, text="Markers to evalate obj functions: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(iden_win,textvariable=evaluate)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(iden_win,text="Enter values",command=msgentry_iden)
	button.grid(row=index,column=0)
	button = Button(iden_win,text="Back",command=iden_win.destroy)
	button.grid(row=index,column=1)

index = 1
label1 = Label(root, text="Problem Type: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,prob,*PROB_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="Restart Solution: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,rest,*YESNO_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="System of Measurements: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,si,*SI_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label = Label(root, text="Regime Type: ")
label.grid(row=index,column=0)
button = Button(root,text="Compressible",command=compressible_window)
button.grid(row=index,column=1)
button = Button(root,text="Incompressible",command=incompressible_window)
button.grid(row=index,column=2)
index+=1

label = Label(root, text="Reference Values: ")
label.grid(row=index,column=0)
button = Button(root,text="Define reference values",command=define_refvalues)
button.grid(row=index,column=1)
index+=1

label = Label(root, text="Timestep attributes: ")
label.grid(row=index,column=0)
button = Button(root,text="Define steady/unsteady parameters",command=define_unsteady)
button.grid(row=index,column=1)
index+=1

label = Label(root, text="Marker attributes: ")
label.grid(row=index,column=0)
button = Button(root,text="Define markers",command=define_markers)
button.grid(row=index,column=1)
button = Button(root,text="Define surface identification ",command=define_iden_markers)
button.grid(row=index,column=2)
index+=1

button = Button(root,text="Write values",command=msgentry_basic)
button.grid(row=index,column=1)
button = Button(root,text="Exit",command=root.quit)
button.grid(row=index,column=3)
index+=1

root.mainloop()
