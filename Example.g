Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");

# load group from CARAT table:
gens := Carat(6, 2377, 10);

# construct a coflasque resolution:
pr := CoflasqueCover(gens);
cr := CoflasqueResolution(pr, gens);

# compute the group Phi(G, M):
result := ComputePhi(gens, cr: num_threads:=1);

Print("\n");
Print(result);
