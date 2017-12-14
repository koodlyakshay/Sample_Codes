import sys
from Tkinter import *

root = Tk()
root.title('SU2 config file')

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

CONVEC_OPTIONS = [
"JST", 
"LAX-FRIEDRICH",
"CUSP", 
"ROE",
"AUSM",
"HLLC",
"TURKEL_PREC",
"MSW"
]

CONVECTURB_OPTIONS = [
"SCALAR UPWIND"
]

TIME_OPTIONS = [
"EULER_IMPLICIT", 
"EULER_EXPLICIT",
"RUNGE_KUTTA"
]

TIMETURB_OPTIONS = [
"EULER_IMPLICIT"
]

SLOPE_OPTIONS = [
"NONE", 
"VENKATAKRISHNAN",
"VENKATAKRISHNAN_WANG", 
"BARTH_JESPERSEN",
"VAN_ALBADA_EDGE"
]

LINSOLVE_OPTIONS = [
"FGMRES", 
"BCGSTAB",
"SMOOTHER_JACOBI", 
"SMOOTHER_ILU",
"SMOOTHER_LUSGS",
"SMOOTHER_LINELET"
]

LINPREC_OPTIONS = [
"ILU",
"JACOBI", 
"LU_SGS",
"LINELET"
]

GRAD_OPTIONS = [
"GREEN_GAUSS",
"WEIGHTED_LEAST_SQUARES"
]

MGCYCLE_OPTIONS = [
"V_CYCLE",
"W_CYCLE", 
"FULLMG_CYCLE"
]

MGLEVEL_OPTIONS = [
"0",
"1", 
"2",
"3"
]

CONVERGE_OPTIONS = [
"CAUCHY",
"RESIDUAL"
]

CAUCHYFUNC_OPTIONS = [
"DRAG",
"LIFT",
"NEARFIELD_PRESS"
]

FORMAT_OPTIONS = [
"SU2",
"CGNS"
]

OPFORMAT_OPTIONS = [
"TECPLOT",
"TECPLOT_BINARY",
"PARAVIEW",
"FIELDVIEW",
"FIELDVIEW_BINARY",
]

GRIDMOV_OPTIONS = [
"RIGID_MOTION",
"DEFORMING",
"MOVING_WALL",
"STEADY_TRANSLATION",
"GUST",
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

convec = StringVar()
convec.set(CONVEC_OPTIONS[0])

convecturb = StringVar()
convecturb.set(CONVECTURB_OPTIONS[0])

slope = StringVar()
slope.set(SLOPE_OPTIONS[0])

linsolve = StringVar()
linsolve.set(LINSOLVE_OPTIONS[0])

linprec = StringVar()
linprec.set(LINPREC_OPTIONS[0])

mg = StringVar()
mg.set(MGCYCLE_OPTIONS[0])

mglevel = StringVar()
mglevel.set(MGLEVEL_OPTIONS[0])

grad = StringVar()
grad.set(GRAD_OPTIONS[0])

cfladapt = StringVar()
cfladapt.set(YESNO_OPTIONS[0])

timefl = StringVar()
timefl.set(TIME_OPTIONS[0])

timetr = StringVar()
timetr.set(TIMETURB_OPTIONS[0])

conv_crit = StringVar()
conv_crit.set(CONVERGE_OPTIONS[0])

cauchy_func = StringVar()
cauchy_func.set(CAUCHYFUNC_OPTIONS[0])

mshfrmt = StringVar()
mshfrmt.set(FORMAT_OPTIONS[0])

opfrmt = StringVar()
opfrmt.set(OPFORMAT_OPTIONS[0])

gridmov = StringVar()
gridmov.set(YESNO_OPTIONS[0])

gridkind = StringVar()
gridkind.set(GRIDMOV_OPTIONS[0])

mshname = StringVar()
mshoutname = StringVar()
mshoutname.set("mesh_out.su2")
solflname = StringVar()
solflname.set("solution_flow.dat")
ophistfile = StringVar()
ophistfile.set("history")
brkfile = StringVar()
brkfile.set("forces_breakdown.dat")
restflname = StringVar()
restflname.set("restart_flow.dat")
restadname = StringVar()
restadname.set("restart_adj.dat")
opfile = StringVar()
opfile.set("flow.dat")
surffile = StringVar()
surffile.set("surface_flow.dat")


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
cfl = StringVar()
cfl.set(10.0)
cflup = StringVar()
cflup.set(0.5)
cfldown = StringVar()
cfldown.set(1.5)
cflmin = StringVar()
cflmin.set(1.25)
cflmax = StringVar()
cflmax.set(50.0)
maxdt = StringVar()
maxdt.set(1E6)
mgpro = StringVar()
mgpro.set(0.75)
mgrest = StringVar()
mgrest.set(0.75)
mgpre= StringVar()
mgpre.set("1,1,1,1")
mgpost= StringVar()
mgpost.set("1,1,1,1")
mgcorrec= StringVar()
mgcorrec.set("1,1,1,1")
enfix = StringVar()
enfix.set(0.0)
relaxfl = StringVar()
relaxfl.set(0.95)
relaxtr = StringVar()
relaxtr.set(0.95)
cflred = StringVar()
cflred.set(1.0)
lfillin = StringVar()
lfillin.set(0)
lerror = StringVar()
lerror.set(1E-6)
liter = StringVar()
liter.set(10)
extiter = StringVar()
extiter.set(99999)
res_red = StringVar()
res_red.set(5)
res_min = StringVar()
res_min.set(-8)
start_iter = StringVar()
start_iter.set(10)
cauchy_elem = StringVar()
cauchy_elem.set(100)
cauchy_eps = StringVar()
cauchy_eps.set(1E-10)
wrtcon = StringVar()
wrtcon.set(1)
wrtcondt = StringVar()
wrtcondt.set(10)
wrtsol = StringVar()
wrtsol.set(1000)
wrtsoldt = StringVar()
wrtsoldt.set(1)

morgx = StringVar()
morgx.set(1.0)
morgy = StringVar()
morgy.set(0.0)
morgz = StringVar()
morgz.set(0.0)

rrx = StringVar()
rrx.set(1.0)
rry = StringVar()
rry.set(0.0)
rrz = StringVar()
rrz.set(0.0)

pampx = StringVar()
pampx.set(1.0)
pampy = StringVar()
pampy.set(0.0)
pampz = StringVar()
pampz.set(0.0)

pomgx = StringVar()
pomgx.set(1.0)
pomgy = StringVar()
pomgy.set(0.0)
pomgz = StringVar()
pomgz.set(0.0)

pphx = StringVar()
pphx.set(1.0)
pphy = StringVar()
pphy.set(0.0)
pphz = StringVar()
pphz.set(0.0)

trrx = StringVar()
trrx.set(1.0)
trry = StringVar()
trry.set(0.0)
trrz = StringVar()
trrz.set(0.0)

plomgx = StringVar()
plomgx.set(1.0)
plomgy = StringVar()
plomgy.set(0.0)
plomgz = StringVar()
plomgz.set(0.0)

plampx = StringVar()
plampx.set(1.0)
plampy = StringVar()
plampy.set(0.0)
plampz = StringVar()
plampz.set(0.0)

machmotion = StringVar()
machmotion.set(0.8)

markmov = StringVar()
markmov.set(NONE)

movmotorg = StringVar()
movmotorg.set(0)

def msgentry_comp():
	f.write("%--------Compressible flow parameters-----------%\n")
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
	f.write("%--------Fluid model parameters-----------%\n")
	f.write("FLUID_MODEL= STANDARD_AIR\n")
	f.write("GAMMA= 1.4\n")
	f.write("CRITICAL_TEMPERATURE= 131.0\n")
	f.write("GAS_CONSTANT= 287.058\n")
	f.write("CRITICAL_PRESSURE= 3588550.0\n")
	f.write("ACCENTRIC_FACTOR= 0.035\n")


def msgentry_flmodel():
	f.write("%--------Fluid Model parameters-----------%\n")
	f.write("FLUID_MODEL= "+flopt.get()+"\n")
	f.write("GAMMA= "+gamma.get()+"\n")
	f.write("CRITICAL_TEMPERATURE= "+t_c.get()+"\n")
	f.write("GAS_CONSTANT= "+gas_cons.get()+"\n")
	f.write("CRITICAL_PRESSURE= "+p_c.get()+"\n")
	f.write("ACCENTRIC_FACTOR= "+a_f.get()+"\n")
	
def msgentry_viscmodel():
	f.write("%--------Viscosity model parameters-----------%\n")
	f.write("VISCOSITY_MODEL= SUTHERLAND\n")
	f.write("MU_REF= "+mu_ref.get()+"\n")
	f.write("MU_T_REF= "+mu_tref.get()+"\n")
	f.write("SUTHERLAND_CONSTANT= "+sut_cons.get()+"\n")
	
def msgentry_consvisc():
	f.write("%--------Viscosity model parameters-----------%\n")
	f.write("VISCOSITY_MODEL= CONSTANT_VISCOSITY\n")
	f.write("MU_CONS= "+mu_cons.get()+"\n")



def msgentry_incomp():
	f.write("%--------Incompressible flow parameters-----------%\n")
	f.write("REGIME_TYPE= INCOMPRESSIBLE\n")
	f.write("FREESTREAM_DENSITY= "+denfs.get()+"\n")
	f.write("FREESTREAM_VELOCITY= ("+velx.get()+","+vely.get()+","+velz.get()+")\n")
	f.write("FREESTREAM_VISCOSITY= "+visc.get()+"\n")

def msgentry_basic():
	f.write("%--------Problem definition-----------%\n")
	f.write("PHYSICAL_PROBLEM= "+prob.get()+"\n")
	f.write("MATH_PROBLEM= DIRECT \n")
	f.write("KIND_TURB_MODEL= "+turb.get()+"\n")
	f.write("RESTART_SOL= "+rest.get()+"\n")
	f.write("SYSTEM_MEASUREMENTS= "+si.get()+"\n")
	
def msgentry_refval():
	f.write("%--------Reference value definition-----------%\n")
	f.write("REF_NONDIMENSIONALISATION= "+refnondim.get()+"\n")
	f.write("REF_LENGTH= "+reflen.get()+"\n")
	f.write("REF_AREA= "+refar.get()+"\n")
	f.write("REF_ORIGIN_MOMENT_X="+refmomx.get()+"\n")
	f.write("REF_ORIGIN_MOMENT_Y="+refmomy.get()+"\n")
	f.write("REF_ORIGIN_MOMENT_Z="+refmomz.get()+"\n")
	f.write("FREESTREAM_VISCOSITY= "+visc.get()+"\n")
	
def msgentry_unsteady():
	f.write("%--------Unsteady simulation parameters-----------%\n")
	f.write("UNSTEADY_SIMULATION= "+unst.get()+"\n")
	f.write("UNST_TIMESTEP= "+unst_delt.get()+"\n")
	f.write("UNST_TIME= "+unst_time.get()+"\n")
	f.write("UNST_CFL_NUMBER= "+unst_cfl.get()+"\n")
	f.write("UNST_INT_ITER= "+n_intiter.get()+"\n")
	f.write("UNST_RESTART_ITER= "+unst_rest.get()+"\n")
	
def msgentry_marker():
	f.write("%--------Boundary conditions-----------%\n")
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
	f.write("%--------Monitioring solution-----------%\n")
	f.write("MARKER_PLOTTING= ("+plotting.get()+")\n")
	f.write("MARKER_MONITORING= ("+monitor.get()+")\n")
	f.write("MARKER_DESIGNING= ("+evaluate.get()+")\n")
	
def msgentry_spacenum():
	f.write("%--------Flow numerical method definition-----------%\n")
	f.write("CONV_NUM_METHOD_FLOW= "+convec.get()+"\n")
	f.write("NUM_METHOD_GRAD= "+grad.get()+"\n")
	f.write("SLOPE_LIMITER_FLOW= "+slope.get()+"\n")
	f.write("TIME_DISCRE_FLOW= "+timefl.get()+"\n")
	f.write("CFL_NUMBER= "+cfl.get()+"\n")
	f.write("CFL_ADAPT= "+cfladapt.get()+"\n")
	f.write("CFL_ADAPT_PARAM= ("+cfldown.get()+","+cflup.get()+","+cflmin.get()+","+cflmax.get()+")\n")
	f.write("MAX_DELTA_TIME= "+maxdt.get()+"\n")
	f.write("ENTROPY_FIX_COEFF= "+enfix.get()+"\n")
	f.write("RELAXATION_FACTOR_FLOW= "+relaxfl.get()+"\n")
	f.write("%--------Turbulent numerical method definition-----------%\n")
	f.write("CONV_NUM_METHOD_TURB= "+convecturb.get()+"\n")
	f.write("TIME_DISCRE_TURB= "+timetr.get()+"\n")
	f.write("RELAXATION_FACTOR_TURB= "+relaxtr.get()+"\n")
	f.write("CFL_REDUCTION_TURB= "+cflred.get()+"\n")
	
	
def msgentry_linsolver():
	f.write("%--------Linear solver parameters-----------%\n")
	f.write("LINEAR_SOLVER= "+linsolve.get()+"\n")
	f.write("LINEAR_SOLVER_PREC= "+linprec.get()+"\n")
	f.write("LINEAR_SOLVER_ILU_FILL_IN= "+lfillin.get()+"\n")
	f.write("LINEAR_SOLVER_ERROR= "+lerror.get()+"\n")
	f.write("LINEAR_SOLVER_ITER= "+liter.get()+"\n")
	
def msgentry_mg():
	f.write("%--------Multigrid parameters-----------%\n")
	f.write("MGLEVEL= "+mglevel.get()+"\n")
	f.write("MGCYCLE= "+mg.get()+"\n")
	f.write("MG_PRE_SMOOTH= ("+mgpre.get()+")\n")
	f.write("MG_POST_SMOOTH= ("+mgpost.get()+")\n")
	f.write("MG_CORRECTION_SMOOTH= ("+mgcorrec.get()+")\n")
	f.write("MG_DAMP_RESTRICTION= "+mgrest.get()+"\n")
	f.write("MG_DAMP_PROLONGATION= "+mgpro.get()+"\n")
	
def msgentry_converge():
	f.write("%--------Convergence criteria-----------%\n")
	f.write("EXTITER= "+extiter.get()+"\n")
	f.write("STARTCONV_ITER= "+start_iter.get()+"\n")


def msgentry_grdnone():
	f.write("%--------Grid Movement-----------%\n")
	f.write("GRID_MOVEMENT= "+gridmov.get()+"\n")
	
	
def msgentry_grdmov():
	f.write("%--------Grid Movement-----------%\n")
	f.write("GRID_MOVEMENT= YES\n")
	f.write("GRID_MOVEMENT_KIND= "+gridkind.get()+"\n")
	f.write("MARKER_MOVING= "+markmov.get()+"\n")
	f.write("MACH_MOTION= "+machmotion.get()+"\n")
	f.write("MOTION_ORIGIN_X="+morgx.get()+"\n")
	f.write("MOTION_ORIGIN_Y="+morgy.get()+"\n")
	f.write("MOTION_ORIGIN_Z="+morgz.get()+"\n")
	f.write("ROTATION_RATE_X="+rrx.get()+"\n")
	f.write("ROTATION_RATE_Y="+rry.get()+"\n")
	f.write("ROTATION_RATE_Z="+rrz.get()+"\n")
	f.write("PITCHING_OMEGA_X="+pomgx.get()+"\n")
	f.write("PITCHING_OMEGA_Y="+pomgy.get()+"\n")
	f.write("PITCHING_OMEGA_Z="+pomgz.get()+"\n")
	f.write("PITCHING_AMPL_X="+pampx.get()+"\n")
	f.write("PITCHING_AMPL_Y="+pampy.get()+"\n")
	f.write("PITCHING_AMPL_Z="+pampz.get()+"\n")
	f.write("PITCHING_PHASE_X="+pphx.get()+"\n")
	f.write("PITCHING_PHASE_Y="+pphy.get()+"\n")
	f.write("PITCHING_PHASE_Z="+pphz.get()+"\n")
	f.write("TRANSLATION_RATE_X="+trrx.get()+"\n")
	f.write("TRANSLATION_RATE_Y="+trry.get()+"\n")
	f.write("TRANSLATION_RATE_Z="+trrz.get()+"\n")
	f.write("PLUNGING_OMEGA_X="+plomgx.get()+"\n")
	f.write("PLUNGING_OMEGA_Y="+plomgy.get()+"\n")
	f.write("PLUNGING_OMEGA_Z="+plomgz.get()+"\n")
	f.write("PLUNGING_AMPL_X="+plampx.get()+"\n")
	f.write("PLUNGING_AMPL_Y="+plampy.get()+"\n")
	f.write("PLUNGING_AMPL_Z="+plampz.get()+"\n")
	f.write("MOVE_MOTION_ORIGIN="+movmotorg.get()+"\n")
	
	


def msgentry_cauchy():
	msgentry_converge()
	f.write("CONV_CRITERIA= CAUCHY\n")
	f.write("CAUCHY_ELEMS= "+cauchy_elem.get()+"\n")
	f.write("CAUCHY_EPS= "+cauchy_eps.get()+"\n")
	f.write("CAUCHY_FUNC_FLOW= "+cauchy_func.get()+"\n")
	
def msgentry_residual():
	msgentry_converge()
	f.write("CONV_CRITERIA= RESIDUAL\n")
	f.write("RESIDUAL_REDUCTION= "+cauchy_elem.get()+"\n")
	f.write("RESIDUAL_MINVAL= "+cauchy_eps.get()+"\n")


def msgentry_ioinfo():
	f.write("%--------Input/Output information-----------%\n")
	f.write("MESH_FILENAME= "+mshname.get()+"\n")
	f.write("MESH_FORMAT= "+mshfrmt.get()+"\n")
	f.write("SOLUTION_FLOW_FILENAME= "+solflname.get()+"\n")
	f.write("OUTPUT_FORMAT= "+opfrmt.get()+"\n")
	f.write("CONV_FILENAME= "+ophistfile.get()+"\n")
	f.write("BREAKDOWN_FILENAME= "+brkfile.get()+"\n")
	f.write("RESTART_FLOW_FILENAME= "+restflname.get()+"\n")
	f.write("VOLUME_FLOW_FILENAME= "+opfile.get()+"\n")
	f.write("SURFACE_FLOW_FILENAME= "+surffile.get()+"\n")
	f.write("WRT_SOL_FREQ= "+wrtsol.get()+"\n")
	f.write("WRT_SOL_FREQ_DUALTIME= "+wrtsoldt.get()+"\n")
	f.write("WRT_CONV_FREQ= "+wrtcon.get()+"\n")
	f.write("WRT_CONV_FREQ_DUALTIME= "+wrtcondt.get()+"\n")
	f.write("WRT_RESIDUALS= NO\n")
	f.write("WRT_LIMITERS= NO\n")
	f.write("WRT_SURFACE= NO\n")
	f.write("LOW_MEMORY_OUTPUT= NO\n")
	f.write("CONSOLE_OUTPUT_VERBOSITY= HIGH\n")
	f.write("WRT_BINARY_RESTART= YES\n")
	f.write("READ_BINARY_RESTART= NO\n")
	
	


f=open(sys.argv[1],"w")



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
	
	button = Button(incomp_win,text="Update cfg file",command=msgentry_incomp)
	button.grid(row=index,column=1)
	index+=1
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
	button = Button(second_win,text="Update cfg file",command=msgentry_comp)
	button.grid(row=index,column=1)
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
	button = Button(consvisc_win,text="Update cfg file",command=msgentry_consvisc)
	button.grid(row=index,column=1)
	index=index+1
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
	
	button = Button(viscmodel_win,text="Update cfg file",command=msgentry_viscmodel)
	button.grid(row=index,column=1)
	index=index+1
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
	
	button = Button(flmodel_win,text="Update cfg file",command=msgentry_flmodel)
	button.grid(row=index,column=1)
	index=index+1
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
	
	button = Button(ref_win,text="Update cfg file",command=msgentry_refval)
	button.grid(row=index,column=1)
	index=index+1
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
	
	button = Button(unsteady_win,text="Update cfg file",command=msgentry_unsteady)
	button.grid(row=index,column=1)
	index=index+1
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
	
	button = Button(marker_win,text="Update cfg file",command=msgentry_marker)
	button.grid(row=index,column=1)
	index=index+1
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
	
	button = Button(iden_win,text="Update cfg file",command=msgentry_iden)
	button.grid(row=index,column=1)
	index=index+1
	button = Button(iden_win,text="Back",command=iden_win.destroy)
	button.grid(row=index,column=1)

def define_spacenum():
	space_win = Toplevel()
	space_win.title("Space Discretization")
	
	index = 1
	
	label1 = Label(space_win, text="Convective numerical method: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(space_win,convec,*CONVEC_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label1 = Label(space_win, text="Gradient method: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(space_win,grad,*GRAD_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label1 = Label(space_win, text="Time discretization: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(space_win,timefl,*TIME_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(space_win, text="CFL number: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(space_win,textvariable=cfl)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(space_win, text="External iterations: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(space_win,textvariable=extiter)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(space_win, text="Start Convergence criteria in iteration: ")
	label2.grid(row=index,columnspan=2)
	entry2 = Entry(space_win,textvariable=start_iter)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label1 = Label(space_win, text="CFL Adapt: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(space_win,cfladapt,*YESNO_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label3 = Label(space_win, text="CFL Adapt parameters")
	label3.grid(row=index,column=0)
	index = index+1
	label3 = Label(space_win, text="Factor down")
	label3.grid(row=index,column=0)
	label3 = Label(space_win, text="Factor up")
	label3.grid(row=index,column=1)
	label3 = Label(space_win, text="CFL min")
	label3.grid(row=index,column=2)
	label3 = Label(space_win, text="CFL max")
	label3.grid(row=index,column=3)
	index = index+1

	entry2 = Entry(space_win,textvariable=cfldown)
	entry2.grid(row=index,column=0)
	entry2 = Entry(space_win,textvariable=cflup)
	entry2.grid(row=index,column=1)
	entry2 = Entry(space_win,textvariable=cflmin)
	entry2.grid(row=index,column=2)
	entry2 = Entry(space_win,textvariable=cflmax)
	entry2.grid(row=index,column=3)
	index=index+1
	
	label2 = Label(space_win, text="Maximum time in local time stepping: ")
	label2.grid(row=index,columnspan=2)
	entry2 = Entry(space_win,textvariable=maxdt)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label2 = Label(space_win, text="Entropy fix coefficient: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(space_win,textvariable=enfix)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(space_win, text="Relaxation coefficient (flow): ")
	label2.grid(row=index,column=0)
	entry2 = Entry(space_win,textvariable=relaxfl)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label1 = Label(space_win, text="Turbulent convective method: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(space_win,convecturb,*CONVECTURB_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label1 = Label(space_win, text="Time discretization(turb): ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(space_win,timetr,*TIMETURB_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(space_win, text="Relaxation coefficient (turb): ")
	label2.grid(row=index,column=0)
	entry2 = Entry(space_win,textvariable=relaxtr)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(space_win, text="CFL Reduction factor (turb): ")
	label2.grid(row=index,column=0)
	entry2 = Entry(space_win,textvariable=cflred)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(space_win,text="Update cfg file",command=msgentry_spacenum)
	button.grid(row=index,column=1)
	index=index+2
	
	label1 = Label(space_win, text="Multigrid parameters: ")
	label1.grid(row=index,column=0)
	button = Button(space_win,text="Define",command=define_mg)
	button.grid(row=index,column=1)
	index+=1
	
	button = Button(space_win,text="Back",command=space_win.destroy)
	button.grid(row=index,column=1)
	
def define_mg():
	mg_win = Toplevel()
	mg_win.title("Multigrid parameters")
	
	index = 1
	label1 = Label(mg_win, text="MG LEVEL: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(mg_win,mglevel,*MGLEVEL_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label1 = Label(mg_win, text="MG CYCLE: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(mg_win,mg,*MGCYCLE_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(mg_win, text="Pre-Smoothing: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(mg_win,textvariable=mgpre)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(mg_win, text="Post-Smoothing: ")
	label2.grid(row=index,column=0)
	entry3 = Entry(mg_win,textvariable=mgpost)
	entry3.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(mg_win, text="Correction smoothing: ")
	label2.grid(row=index,column=0)
	entry4 = Entry(mg_win,textvariable=mgcorrec)
	entry4.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(mg_win, text="Restriction damping factor: ")
	label2.grid(row=index,column=0)
	entry4 = Entry(mg_win,textvariable=mgrest)
	entry4.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(mg_win, text="Prolongation damping factor: ")
	label2.grid(row=index,column=0)
	entry4 = Entry(mg_win,textvariable=mgpro)
	entry4.grid(row=index,column=1)
	index=index+1
	
	button = Button(mg_win,text="Update cfg file",command=msgentry_mg)
	button.grid(row=index,column=1)
	index=index+1
	button = Button(mg_win,text="Back",command=mg_win.destroy)
	button.grid(row=index,column=1)
	
def define_linsol():
	lsol_win=Toplevel()
	lsol_win.title("Linear solver definition")
	index = 1
	
	label1 = Label(lsol_win, text="Linear solver: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(lsol_win,linsolve,*LINSOLVE_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label1 = Label(lsol_win, text="Preconditioner: ")
	label1.grid(row=index,column=0)
	opt1 = OptionMenu(lsol_win,linprec,*LINPREC_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(lsol_win, text="ILU Preconditioner fill-in: ")
	label2.grid(row=index,columnspan=2)
	entry2 = Entry(lsol_win,textvariable=lfillin)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label2 = Label(lsol_win, text="Minimum error for implicit: ")
	label2.grid(row=index,columnspan=2)
	entry2 = Entry(lsol_win,textvariable=lerror)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label2 = Label(lsol_win, text="Max iterations of linear solver for implicit: ")
	label2.grid(row=index,columnspan=2)
	entry2 = Entry(lsol_win,textvariable=liter)
	entry2.grid(row=index,column=2)
	index=index+1
	
	button = Button(lsol_win,text="Update cfg file",command=msgentry_linsolver)
	button.grid(row=index,column=1)
	index=index+1
	button = Button(lsol_win,text="Back",command=lsol_win.destroy)
	button.grid(row=index,column=1)
	
def define_cauchy():
	cauchy_win=Toplevel()
	cauchy_win.title("Cauchy convergence parameters")
	index = 1	

	label2 = Label(cauchy_win, text="Cauchy Elements: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(cauchy_win,textvariable=cauchy_elem)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(cauchy_win, text="Epsilon: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(cauchy_win,textvariable=cauchy_eps)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(cauchy_win, text="Cauchy functions: ")
	label2.grid(row=index,column=0)
	opt1 = OptionMenu(cauchy_win,cauchy_func,*CAUCHYFUNC_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	button = Button(cauchy_win,text="Update cfg file",command=msgentry_cauchy)
	button.grid(row=index,column=1)
	index=index+1
	button = Button(cauchy_win,text="Back",command=cauchy_win.destroy)
	button.grid(row=index,column=1)


def define_residual():
	res_win=Toplevel()
	res_win.title("Residual convergence parameters")
	index = 1	
	label2 = Label(res_win, text="Residual reduction: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(res_win,textvariable=res_red)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(res_win, text="Min value of residual: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(res_win,textvariable=res_min)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(res_win,text="Update cfg file",command=msgentry_residual)
	button.grid(row=index,column=1)
	index=index+1
	button = Button(res_win,text="Back",command=res_win.destroy)
	button.grid(row=index,column=1)


def define_ioinfo():
	io_win=Toplevel()
	io_win.title("Input/Output information")
	index = 1
	
	label2 = Label(io_win, text="Mesh input file: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=mshname)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Mesh format: ")
	label2.grid(row=index,column=0)
	opt1 = OptionMenu(io_win,mshfrmt,*FORMAT_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Mesh output file: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=mshoutname)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Restart flow input file: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=solflname)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Output file format: ")
	label2.grid(row=index,column=0)
	opt1 = OptionMenu(io_win,opfrmt,*OPFORMAT_OPTIONS)
	opt1.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Convergence history file(w/o extension): ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=ophistfile)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Forces breakdown file: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=brkfile)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Restart flow output file: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=restflname)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Flow variable output file: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=opfile)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Surface flow coefficient output file: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=surffile)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Writing solution file frequency: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=wrtsol)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Solution file frequency(dual time): ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=wrtsoldt)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Convergency history frequency: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=wrtcon)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label2 = Label(io_win, text="Convergency history frequency(dual time): ")
	label2.grid(row=index,column=0)
	entry2 = Entry(io_win,textvariable=wrtcondt)
	entry2.grid(row=index,column=1)
	index=index+1
	
	button = Button(io_win,text="Update cfg file",command=msgentry_ioinfo)
	button.grid(row=index,column=1)
	index=index+1
	button = Button(io_win,text="Back",command=io_win.destroy)
	button.grid(row=index,column=1)



def define_grd():
	grd_win=Toplevel()
	grd_win.title("Input/Output information")
	index = 1
	
	
	label2 = Label(grd_win, text="Mach Motion: ")
	label2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=machmotion)
	entry2.grid(row=index,column=1)
	index=index+1
	
	label3 = Label(grd_win, text="Motion Origin(x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(grd_win,textvariable=morgx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=morgy)
	entry2.grid(row=index,column=1)
	entry2 = Entry(grd_win,textvariable=morgz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label3 = Label(grd_win, text="Rotation Rate(x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(grd_win,textvariable=rrx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=rry)
	entry2.grid(row=index,column=1)
	entry2 = Entry(grd_win,textvariable=rrz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label3 = Label(grd_win, text="Pitching Angular freq(rad/s, x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(grd_win,textvariable=pomgx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=pomgy)
	entry2.grid(row=index,column=1)
	entry2 = Entry(grd_win,textvariable=pomgz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label3 = Label(grd_win, text="Pitching Amplitude, degrees(x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(grd_win,textvariable=pampx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=pampy)
	entry2.grid(row=index,column=1)
	entry2 = Entry(grd_win,textvariable=pampz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label3 = Label(grd_win, text="Pitching phase offset, degrees(x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(grd_win,textvariable=pphx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=pphy)
	entry2.grid(row=index,column=1)
	entry2 = Entry(grd_win,textvariable=pphz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label3 = Label(grd_win, text="Translational velocity, m/s(x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(grd_win,textvariable=trrx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=trry)
	entry2.grid(row=index,column=1)
	entry2 = Entry(grd_win,textvariable=trrz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label3 = Label(grd_win, text="Plunging Angular freq(rad/s, x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(grd_win,textvariable=plomgx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=plomgy)
	entry2.grid(row=index,column=1)
	entry2 = Entry(grd_win,textvariable=plomgz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label3 = Label(grd_win, text="Plunging Amplitude, degrees(x,y,z)")
	label3.grid(row=index,column=0)
	index = index+1

	entry2 = Entry(grd_win,textvariable=plampx)
	entry2.grid(row=index,column=0)
	entry2 = Entry(grd_win,textvariable=plampy)
	entry2.grid(row=index,column=1)
	entry2 = Entry(grd_win,textvariable=plampz)
	entry2.grid(row=index,column=2)
	index=index+1
	
	label2 = Label(grd_win, text="Move motion origin for moving marker (1 or 0): ")
	label2.grid(row=index,columnspan=2)
	entry2 = Entry(grd_win,textvariable=movmotorg)
	entry2.grid(row=index,column=2)
	index=index+1
	
	button = Button(grd_win,text="Update cfg file",command=msgentry_grdmov)
	button.grid(row=index,column=1)
	index=index+1
	button = Button(grd_win,text="Back",command=grd_win.destroy)
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

label1 = Label(root, text="Turbulence Model: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,turb,*TURB_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

label1 = Label(root, text="System of Measurements: ")
label1.grid(row=index,column=0)
opt1 = OptionMenu(root,si,*SI_OPTIONS)
opt1.grid(row=index,column=1)
index=index+1

button = Button(root,text="Update cfg file",command=msgentry_basic)
button.grid(row=index,column=1)
index=index+1

label = Label(root, text="Regime Type: ")
label.grid(row=index,column=0)
button = Button(root,text="Compressible",command=compressible_window)
button.grid(row=index,column=1)
button = Button(root,text="Incompressible",command=incompressible_window)
button.grid(row=index,column=2)
index+=1

label = Label(root, text="Grid Movement: ")
label.grid(row=index,column=0)
button = Button(root,text="None",command=msgentry_grdnone)
button.grid(row=index,column=1)
button = Button(root,text="Define",command=define_grd)
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

label = Label(root, text="Numerical Methods: ")
label.grid(row=index,column=0)
button = Button(root,text="Discretization",command=define_spacenum)
button.grid(row=index,column=1)
button = Button(root,text="Linear solver",command=define_linsol)
button.grid(row=index,column=2)
index+=1

label = Label(root, text="Marker attributes: ")
label.grid(row=index,column=0)
button = Button(root,text="Define markers",command=define_markers)
button.grid(row=index,column=1)
button = Button(root,text="Define surface identification ",command=define_iden_markers)
button.grid(row=index,column=2)
index+=1

label = Label(root, text="Convergence parameters: ")
label.grid(row=index,column=0)
button = Button(root,text="Cauchy",command=define_cauchy)
button.grid(row=index,column=1)
button = Button(root,text="Residual ",command=define_residual)
button.grid(row=index,column=2)
index+=1

label = Label(root, text="Input/Output information: ")
label.grid(row=index,column=0)
button = Button(root,text="Define",command=define_ioinfo)
button.grid(row=index,column=1)

index+=1

button = Button(root,text="Exit",command=root.quit)
button.grid(row=index,column=2)
index+=1

root.mainloop()
