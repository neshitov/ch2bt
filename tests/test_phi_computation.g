#######################################################
# Tests for computation of Phi constructions
#######################################################

Read("./ComputePhi.g");
Read("./CaratReader.g");
Read("./CoflasqueCover.g");

u := [1, 2, 3, 4];
v := [0, 1, 3, -2];

# exterior product in basis e1e2 e1e3 e1e4 e2e3 e2e4 e3e4
product := [1, 3, -2, 3, -8, -18];

if  ExteriorProduct(u, v) <> product then
  Error("Exterior product error");
fi;

Print("pass\n");

m := [[1,2,3],
      [4,0,1],
      [5,1,2]];

# exterior square in basis e1e2 e1e3 e1e4 e2e3 e2e4 e3e4
m_ex_sq := [[-8, -11, 2],
            [-9, -13, 1],
            [4, 3, -1]];

if  ExteriorSquare([m]) <> [m_ex_sq] then
  Print(ExteriorSquare([m]));
  Error("Exterior square error");
fi;

Print("pass\n");

gens := Carat(5, 697, 4);
cr := CoflasqueResolution(CoflasqueCover(gens), gens);
mapk := MapK(cr, Length(gens));
l2_action := ExteriorSquare(cr.actionC);
for i in [1 .. Length(gens)] do
  if mapk * (Z(2) * l2_action[i])  <> (Z(2) * gens[i]) * mapk then
    Error("Map k not equivariant");
  fi;
od;

Print("pass\n");

gens := Carat(6, 2377, 10);
cr := CoflasqueResolution(CoflasqueCover(gens), gens);
result := ComputePhi(gens, cr);
if result <> rec( HM2_rank := 3, Phi_rank := 1, im_HL_rank := 1, im_HP2_rank := 0 ) then
  Error("Phi computation error");
fi;

Print("\npass\n");
