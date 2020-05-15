


# default parameters to calculate in primer3.
mv_conc=50
dv_conc=8
dntp_conc=1.4

# distance parameters
# ====F3=g1=F2=g2=F1=g3=B1=g4=B2=g5=B3====
#           =========g6=========
g1 = (3,60)
g2 = (25,40)
g3 = (3,40)
g4 = (25,40)
g5 = (3,60)
g6 = (120,160)

# primer length parameters
P3L = (19,21)
P2L = (19,21)
P1L = (20,22)
LPL = (18,25)

# Primer Tm parameters
P3Tm = (62,63)
P2Tm = (62,63)
P1Tm = (67,68)
LPTm = (65.5,66.5)
# Primer GC ratio
GCratio = (0.4,0.65)

# inclusitivity of core Fragments
PInclu = 0.99

# inclusivity of Loop primers
LoopInclu = 0.99

# 3'End stability , simplify as 2 out 5 n.t. must be C or G at 3' end.
E3 = 2

# 3'End self complement length threshod.
ESC = 4

# primer dimer check Tm thereshold
PrimerDimerTm = 25

# primer hairpin structure dG threshold
HairpindG = -2.5
# internal Loop hairpin
LoopHairpindG = -2.5

# bat coV crossreactivity
BatHomology = 2
BatHomologyEnd = 3

# Adjust steps if anything didn't find.
AdjustStep = 3
# primer set with the same F2 threshold
F2_CountThreshold = 50
F1_CountThreshold = 25
B1c_CountThreshold = 15
B2c_CountThreshold = 8
B3c_CountThreshold = 4
