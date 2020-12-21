# The following code is used in
# "Rationality Problem for Algebraic Tori"
# by Akinari Hoshi and Aiichi Yamasaki
# Written by Aiichi Yamasaki
# The files in crystdat.zip are required.

carat_folder := "./carat_tables";
files := [
  Concatenation(carat_folder, "/cryst1.txt"),
  Concatenation(carat_folder, "/cryst2.txt"),
  Concatenation(carat_folder, "/cryst3.txt"),
  Concatenation(carat_folder, "/cryst4.txt"),
  Concatenation(carat_folder, "/cryst5.txt"),
  Concatenation(carat_folder, "/cryst6.txt"),
  Concatenation(carat_folder, "/caratchpol.txt"),
];

for file in files do
  Read(file);
od;


cryst:=[cryst1,cryst2,cryst3,cryst4,cryst5,cryst6];
chpol:=[chpol1,chpol2,chpol3,chpol4,chpol5,chpol6];

CharacteristicCyclotomicPolynomial:= function(m)
    local d,l,i;
    d:=Length(m);
    l:=[];
    for i in DivisorsInt(Order(m)) do
        l:=Concatenation(l,List([RankMat(m-E(i)*IdentityMat(d))+1..d],x->i));
    od;
    return l;
end;

CharPolyList:= function(g)
    local c;
    c:=List(ConjugacyClasses(g),Representative);
    return SortedList(List(c,CharacteristicCyclotomicPolynomial));
end;

CharPolySubgroupsList:= function(g)
    local h;
    h:=List(ConjugacyClassesSubgroups(g),Representative);
    return SortedList(List(h,CharPolyList));
end;

CaratQClass:= function(G)
    local d,ch,q;
    d:=Length(Identity(G));
    ch:=CharPolyList(G);
    q:=Position(chpol[d],ch);
    if d=6 and (q in chpol6dup) then
        ch:=CharPolySubgroupsList(G);
        q:=chpol6dup[Position(chpol6sub,ch)];
    fi;
    return [d,q];
end;
