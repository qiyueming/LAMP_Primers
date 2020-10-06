# default parameters to calculate in primer3.
mv_conc=50
dv_conc=8
dntp_conc=1.4

# distance parameters
# ====F3=g1=F2=g2=F1=g3=B1=g4=B2=g5=B3====
#           =========g6=========
g1 = (1,60)
g2 = (20,40)
g3 = (0,40)
g4 = (20,40)
g5 = (1,60)
g6 = (120,160)

# primer length parameters
P3L = (18,25)
P2L = (18,25)
P1L = (18,25)
LPL = (18,25)

# Primer Tm parameters
P3Tm = (61.5,63.5)
P2Tm = (61.5,63.5)
P1Tm = (66.5,68.5)
LPTm = (65,67)
# Primer GC ratio
GCratio = (0.4,0.65)

# inclusitivity of core Fragments
PInclu = 0.999

# inclusivity of Loop primers
LoopInclu = 0.999

# 3'End stability , simplify as 2 out 5 n.t. must be C or G at 3' end.
E3 = 2

# 3'End self complement length threshod.
ESC = 4

# primer dimer check Tm thereshold
PrimerDimerTm = 30

# primer hairpin structure dG threshold
HairpindG = -4
# internal Loop hairpin
LoopHairpindG = -4

# bat coV crossreactivity
BatHomology = 1
# how many bases at end that need mutaiton.
BatHomologyEnd = 0
# how many fragments in primer set need to satisfy homology.
BatHomologyCount = 0

# Adjust steps if anything didn't find.
AdjustStep = 1
# primer set with the same F2 threshold
F2_CountThreshold = 50
F1_CountThreshold = 25
B1c_CountThreshold = 20
B2c_CountThreshold = 8
B3c_CountThreshold = 4
# max number of the ssame fragment.
SAME_Fragment_Threshold = 10 
PRIMER_DESIGN_METHOD = "main_Inclusivity"
