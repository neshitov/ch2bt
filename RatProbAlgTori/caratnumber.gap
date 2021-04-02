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
  Concatenation(carat_folder, "/H1cryst.txt"),
];



for file in files do
  Read(file);
od;


cryst:=[cryst1,cryst2,cryst3,cryst4,cryst5,cryst6];
chpol:=[chpol1,chpol2,chpol3,chpol4,chpol5,chpol6];
h1glnz:=[h1gl1z,h1gl2z,h1gl3z,h1gl4z,h1gl5z,h1gl6z];

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

Hminus1:= function(g)
    local m,gg,i,s,r;
    m:=[];
    gg:=GeneratorsOfGroup(g);
    if gg=[] then
        return [];
    else
        for i in gg do
            m:=Concatenation(m,i-Identity(g));
        od;
        s:=SmithNormalFormIntegerMat(m);
        r:=Rank(s);
        return List([1..r],x->s[x][x]);
    fi;
end;

H0:= function(g)
    local m,s,r;
    m:=Sum(g);
    s:=SmithNormalFormIntegerMat(m);
    r:=Rank(s);
    return List([1..r],x->s[x][x]);
end;

H1:= function(g)
    local m,gg,i,s,r;
    m:=[];
    gg:=GeneratorsOfGroup(g);
    if gg=[] then
        return [];
    else
        for i in gg do
            m:=Concatenation(m,TransposedMat(i)-Identity(g));
        od;
        m:=TransposedMat(m);
        s:=SmithNormalFormIntegerMat(m);
        r:=Rank(s);
        return List([1..r],x->s[x][x]);
    fi;
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


CaratZClass:= function(G)
    local cat,z,zz,h,hh,zg;
    cat:=CaratQClass(G);
    hh:=h1glnz[cat[1]][cat[2]];
    if Length(hh)=1 then
        return Concatenation(cat,[1]);
    else
        h:=[Hminus1(G),H0(G),H1(G)];
        zz:=Filtered([1..Length(hh)],x->hh[x]=h);
        if Length(zz)=1 then
            return Concatenation(cat,zz);
        else
            for z in zz do
                zg:=Group(cryst[cat[1]][cat[2]][z],IdentityMat(cat[1]));
                if RepresentativeAction(GL(cat[1],Integers),G,zg) in GL(cat[1],Integers) then
                    return Concatenation(cat,[z]);
                fi;
            od;
        fi;
    fi;
end;
