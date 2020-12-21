# The following code is used in
# "Rationality Problem for Algebraic Tori"
# by Akinari Hoshi and Aiichi Yamasaki
# Written by Aiichi Yamasaki

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

Z0lattice:= function(g)
    local gg,m,i;
    gg:=GeneratorsOfGroup(g);
    if gg=[] then
        return Identity(g);
    else
        m:=[];
        for i in gg do
            m:=Concatenation(m,TransposedMat(i)-Identity(g));
        od;
        m:=TransposedMat(m);
        return NullspaceIntMat(m);
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

ConjugacyClassesSubgroups2:= function(g)
    Reset(GlobalMersenneTwister);
    Reset(GlobalRandomSource);
    return ConjugacyClassesSubgroups(g);
end;

ConjugacyClassesSubgroupsFromPerm:= function(g)
    local iso,h,i;
    Reset(GlobalMersenneTwister);
    Reset(GlobalRandomSource);
    iso:=IsomorphismPermGroup(g);
    h:=ConjugacyClassesSubgroups2(Range(iso));
    h:=List(h,Representative);
    h:=List(h,x->PreImage(iso,x));
    return h;
end;

IsFlabby:= function(g)
    local h;
    h:=List(ConjugacyClassesSubgroups2(g),Representative);
    return ForAll(h,x->Product(Hminus1(x))=1);
end;

IsCoflabby:= function(g)
    local h;
    h:=List(ConjugacyClassesSubgroups2(g),Representative);
    return ForAll(h,x->Product(H1(x))=1);
end;

ReduceCoflabbyResolutionBase:= function(g,hh,mi)
    local o,oo,z0,mi2;
    oo:=Orbits(g,mi);
    z0:=List(hh,Z0lattice);
    for o in oo do
       mi2:=Filtered(mi,x->(x in o)=fail);
       if ForAll([1..Length(hh)],i->LatticeBasis(List(Orbits(hh[i],mi2),Sum))=z0[i]) then
           mi:=mi2;
       fi;
    od;
    return mi;
end;

FindCoflabbyResolutionBase:= function(g,hh)
    local d,mi,h,z0,ll,i,o;
    d:=Length(Identity(g));
    mi:=[];
    for h in hh do
        z0:=Z0lattice(h);
        ll:=LatticeBasis(List(Orbits(h,mi),Sum));
        for i in z0 do
            if LatticeBasis(Concatenation(ll,[i]))<>ll then
                o:=Orbit(g,i);
                mi:=Concatenation(mi,o);
                o:=List(Orbits(h,o),Sum);
                ll:=LatticeBasis(Concatenation(ll,o));
            fi;
        od;
    od;
    return    ReduceCoflabbyResolutionBase(g,hh,mi);
end;

FlabbyResolution:= function(g)
    local tg,gg,d,th,mi,ms,o,r,gg1,gg2,v1,mg,img;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    o:=IdentityMat(r);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return rec(injection:=TransposedMat(mi),
                   surjection:=NullMat(r,0),
                   actionP:=TransposedMatrixGroup(Group(gg1,o))
        );
    else
        ms:=NullspaceIntMat(mi);
        v1:=NullspaceIntMat(TransposedMat(ms));
        mg:=Concatenation(v1,ms);
        img:=mg^-1;
        gg2:=List(gg1,x->mg*x*img);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        return rec(injection:=TransposedMat(mi),
                   surjection:=TransposedMat(ms),
                   actionP:=TransposedMatrixGroup(Group(gg1)),
                   actionF:=TransposedMatrixGroup(Group(gg2))
        );
    fi;
end;

ModpTest:= function(g,p)
    local tg,gg,d,th,mi,z0,ll,h,r,i,j,k,gg1,g1,oo,iso,l1,l2,m1,m2;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->Permutation(x,mi));
    gg:=GeneratorsOfGroup(g);
    gg1:=List(gg1,x->x^-1);
    g1:=Group(gg1,());
    mi:=TransposedMat(mi);
    iso:=GroupHomomorphismByImagesNC(g1,g,gg1,gg);
    oo:=OrbitsDomain(g1,[1..r]);
    m1:=[];
    for i in oo do
        h:=Stabilizer(g1,i[1]);
        ll:=List(RightCosetsNC(g1,h),Representative);
        l1:=List(ll,x->i[1]^x);
        l2:=List(ll,x->Image(iso,x));
        z0:=Z0lattice(Image(iso,h));
        for j in z0 do
            m2:=NullMat(r,d);
            for k in [1..Length(ll)] do
                m2[l1[k]]:=j*l2[k];
            od;
            Add(m1,Flat(mi*m2));
        od;
    od;
    m2:=Concatenation(m1,[Flat(IdentityMat(d))]);
    m1:=m1*Z(p);
    m2:=m2*Z(p);
    return RankMat(m1)=RankMat(m2);
end;

SplitTest:= function(g)
    local tg,gg,d,th,mi,z0,ll,h,r,i,j,k,gg1,g1,oo,iso,l1,l2,m1,m2;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->Permutation(x,mi));
    gg:=GeneratorsOfGroup(g);
    gg1:=List(gg1,x->x^-1);
    g1:=Group(gg1,());
    mi:=TransposedMat(mi);
    iso:=GroupHomomorphismByImagesNC(g1,g,gg1,gg);
    oo:=OrbitsDomain(g1,[1..r]);
    m1:=[];
    for i in oo do
        h:=Stabilizer(g1,i[1]);
        ll:=List(RightCosetsNC(g1,h),Representative);
        l1:=List(ll,x->i[1]^x);
        l2:=List(ll,x->Image(iso,x));
        z0:=Z0lattice(Image(iso,h));
        for j in z0 do
            m2:=NullMat(r,d);
            for k in [1..Length(ll)] do
                m2[l1[k]]:=j*l2[k];
            od;
            Add(m1,Flat(mi*m2));
        od;
    od;
    m2:=Concatenation(m1,[Flat(IdentityMat(d))]);
    m1:=NullspaceIntMat(m2);
    return Gcd(TransposedMat(m1)[Length(m2)])=1;
end;

IsInvertibleF:= function(g)
    local tg,gg,d,th,mi,mi2,ms,z0,ll,h,r,i,j,k,gg1,gg2,g1,g2,oo,iso,ker,tg2,th2,h2,l1,l2,v1,m1,m2;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v1:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v1,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        for h in th2 do
            if Product(Hminus1(h))>1 then
                return false;
            fi;
        od;
        for h in th do
            if Product(Hminus1(h))>1 then
                d:=r-d;
                mi:=FindCoflabbyResolutionBase(tg2,th2);
                r:=Length(mi);
                gg1:=List(gg2,x->Permutation(x,mi));
                gg2:=List(gg2,TransposedMat);
                g2:=Group(gg2);
                gg1:=List(gg1,x->x^-1);
                g1:=Group(gg1);
                mi:=TransposedMat(mi);
                iso:=GroupHomomorphismByImagesNC(g1,g2,gg1,gg2);
                oo:=OrbitsDomain(g1,[1..r]);
                m1:=[];
                for i in oo do
                    h:=Stabilizer(g1,i[1]);
                    ll:=List(RightCosetsNC(g1,h),Representative);
                    l1:=List(ll,x->i[1]^x);
                    l2:=List(ll,x->Image(iso,x));
                    z0:=Z0lattice(Image(iso,h));
                    for j in z0 do
                        m2:=NullMat(r,d);
                        for k in [1..Length(ll)] do
                            m2[l1[k]]:=j*l2[k];
                        od;
                        Add(m1,Flat(mi*m2));
                    od;
                od;
                m2:=Concatenation(m1,[Flat(IdentityMat(d))]);
                m1:=NullspaceIntMat(m2);
                return Gcd(TransposedMat(m1)[Length(m2)])=1;
            fi;
        od;
        gg1:=List(gg,x->Permutation(x,mi));
        gg:=GeneratorsOfGroup(g);
        gg1:=List(gg1,x->x^-1);
        g1:=Group(gg1);
        mi:=TransposedMat(mi);
        iso:=GroupHomomorphismByImagesNC(g1,g,gg1,gg);
        oo:=OrbitsDomain(g1,[1..r]);
        m1:=[];
        for i in oo do
            h:=Stabilizer(g1,i[1]);
            ll:=List(RightCosetsNC(g1,h),Representative);
            l1:=List(ll,x->i[1]^x);
            l2:=List(ll,x->Image(iso,x));
            z0:=Z0lattice(Image(iso,h));
            for j in z0 do
                m2:=NullMat(r,d);
                for k in [1..Length(ll)] do
                    m2[l1[k]]:=j*l2[k];
                od;
                Add(m1,Flat(mi*m2));
            od;
        od;
        m2:=Concatenation(m1,[Flat(IdentityMat(d))]);
        m1:=NullspaceIntMat(m2);
        return Gcd(TransposedMat(m1)[Length(m2)])=1;
    fi;
end;

flfl:= function(g)
    local tg,gg,d,th,mi,ms,r,gg1,gg2,v1,mg,img,tg2,iso,ker;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return [];
    else
        ms:=NullspaceIntMat(mi);
        v1:=NullspaceIntMat(TransposedMat(ms));
        mg:=Concatenation(v1,ms);
        img:=mg^-1;
        gg2:=List(gg1,x->mg*x*img);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th:=List(Filtered(th,x->IsSubset(x,ker)),x->Image(iso,x));
        tg:=tg2;
        gg:=gg2;
        d:=r-d;
        mi:=FindCoflabbyResolutionBase(tg,th);
        r:=Length(mi);
        gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
        if r=d then
            return [];
        else
            ms:=NullspaceIntMat(mi);
            v1:=NullspaceIntMat(TransposedMat(ms));
            mg:=Concatenation(v1,ms);
            img:=mg^-1;
            gg2:=List(gg1,x->mg*x*img);
            gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
            return TransposedMatrixGroup(Group(gg2));
        fi;
    fi;
end;

PossibilityOfStablyPermutationH:= function(g,hh)
    local gg,hg,hgg,hom,c,h,m,m1,m2,v,h0,og,oh,p,e;
    gg:=GeneratorsOfGroup(g);
    hg:=List(hh,x->RightCosets(g,x));
    hgg:=List(hg,x->List(gg,y->Permutation(y,x,OnRight)));
    hom:=List(hgg,x->GroupHomomorphismByImages(g,Group(x,()),gg,x));
    c:=List(ConjugacyClasses(g),Representative);
    og:=Order(g);
    m:=List(c,x->List([1..Length(hh)],y->og/Order(hh[y])-NrMovedPoints(Image(hom[y],x))));
    v:=List(c,Trace);
    for h in hh do
        h0:=H0(h);
        oh:=Order(h);
        m1:=List([1..Length(hh)],x->List(OrbitLengths(Image(hom[x],h),[1..og/Order(hh[x])]),y->oh/y));
        Add(m,List(m1,Length));
        Add(v,Length(h0));
        if oh>1 then
            for p in Set(FactorsInt(oh)) do
                for e in [1..PadicValuation(oh,p)] do
                    Add(m,List(m1,x->Number(x,y->PadicValuation(y,p)=e)));
                    Add(v,Number(h0,x->PadicValuation(x,p)=e));
                od;
            od;
        fi;
    od;
    m:=TransposedMat(m);
    m:=Concatenation(m,[v]);
    return NullspaceIntMat(m);
end;

PossibilityOfStablyPermutationM:= function(g)
    local hh;
    hh:=List(ConjugacyClassesSubgroups2(g),Representative);
    return PossibilityOfStablyPermutationH(g,hh);
end;

PossibilityOfStablyPermutationF:= function(g)
    local tg,gg,d,th,mi,ms,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return [0,1];
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return PossibilityOfStablyPermutationH(g2,h2);
    fi;
end;

CosetRepresentation:= function(g,h)
    local gg,hg,og;
    gg:=GeneratorsOfGroup(g);
    hg:=SortedList(RightCosets(g,h));
    og:=Length(hg);
    return List(gg,x->PermutationMat(Permutation(x,hg,OnRight),og));
end;

StablyPermutationCheckH:= function(g,hh,c1,c2)
    local gg,g1,g2,dx,m,i,j,d;
    gg:=List(hh,x->CosetRepresentation(g,x));
    Add(gg,GeneratorsOfGroup(g));
    g1:=[];
    g2:=[];
    for i in [1..Length(gg[1])] do
        m:=[];
        for j in [1..Length(gg)] do
            m:=Concatenation(m,List([1..c1[j]],x->gg[j][i]));
        od;
        Add(g1,DirectSumMat(m));
        m:=[];
        for j in [1..Length(gg)] do
            m:=Concatenation(m,List([1..c2[j]],x->gg[j][i]));
        od;
        Add(g2,DirectSumMat(m));
    od;
    d:=Length(g1[1]);
    if d<>Length(g2[1]) then
        return fail;
    else
        return RepresentativeAction(GL(d,Integers),Group(g1),Group(g2));
    fi;
end;

StablyPermutationMCheck:= function(g,c1,c2)
    local h;
    h:=List(ConjugacyClassesSubgroups2(g),Representative);
    return StablyPermutationCheckH(g,h,c1,c2);
end;

StablyPermutationFCheck:= function(g,c1,c2)
    local tg,gg,d,th,mi,ms,o,h,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return StablyPermutationCheckH(g2,h2,c1,c2);
    fi;
end;

Plist:= function(l)
    return List(l,x->Maximum([x,0]));
end;

Nlist:= function(l)
    return List(l,x->Maximum([-x,0]));
end;

TransformationMat:= function(l1,l2)
    local d1,d2,l,m,p,i,j;
    d1:=Length(l1[1]);
    d2:=Length(l2[1]);
    l:=Length(l1);
    m:=[];
    for i in [1..d1] do
        for j in [1..d2] do
            p:=NullMat(d1,d2);
            p[i][j]:=1;
            Add(m,Flat(List([1..l],x->l1[x]*p-p*l2[x])));
        od;
    od;
    p:=NullspaceIntMat(m);
    return List(p,x->List([1..d1],y->x{[(y-1)*d2+1..y*d2]}));
end;

StablyPermutationCheckHP:= function(g,hh,c1,c2)
    local gg,l,m,m1,m2,tm,d1,d2,s1,s2,i,j,k;
    gg:=List(hh,x->CosetRepresentation(g,x));
    Add(gg,GeneratorsOfGroup(g));
    l:=List([1..Length(gg)],x->Length(gg[x][1]));
    d1:=Sum([1..Length(gg)],x->c1[x]*l[x]);
    d2:=Sum([1..Length(gg)],x->c2[x]*l[x]);
    m:=[];
    s1:=0;
    for i in [1..Length(gg)] do
        if c1[i]>0 then
            m1:=[];
            s2:=0;
            for j in [1..Length(gg)] do
                if c2[j]>0 then
                    tm:=TransformationMat(gg[i],gg[j]);
                    for k in [1..c2[j]] do
                        m2:=List(tm,x->TransposedMat(Concatenation(
                          [NullMat(s2,l[i]),TransposedMat(x),
                          NullMat(d2-s2-l[j],l[i])])));
                        m1:=Concatenation(m1,m2);
                        s2:=s2+l[j];
                    od;
                fi;
            od;
            m1:=LatticeBasis(List(m1,Flat));
            m1:=List(m1,x->List([1..l[i]],y->x{[(y-1)*s2+1..y*s2]}));
            for k in [1..c1[i]] do
                m:=Concatenation(m,List(m1,x->Concatenation(
                  [NullMat(s1,d2),x,NullMat(d1-s1-l[i],d2)])));
                s1:=s1+l[i];
            od;
        fi;
    od;
    return m;
end;

StablyPermutationMCheckP:= function(g,c1,c2)
    local h;
    h:=List(ConjugacyClassesSubgroups2(g),Representative);
    return StablyPermutationCheckHP(g,h,c1,c2);
end;

StablyPermutationFCheckP:= function(g,c1,c2)
    local tg,gg,d,th,mi,ms,o,h,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return StablyPermutationCheckHP(g2,h2,c1,c2);
    fi;
end;

StablyPermutationCheckHMat:= function(g,hh,c1,c2,p)
    local gg,g1,g2,dx,m,i,j,d;
    gg:=List(hh,x->CosetRepresentation(g,x));
    Add(gg,GeneratorsOfGroup(g));
    g1:=[];
    g2:=[];
    for i in [1..Length(gg[1])] do
        m:=[];
        for j in [1..Length(gg)] do
            m:=Concatenation(m,List([1..c1[j]],x->gg[j][i]));
        od;
        Add(g1,DirectSumMat(m));
        m:=[];
        for j in [1..Length(gg)] do
            m:=Concatenation(m,List([1..c2[j]],x->gg[j][i]));
        od;
        Add(g2,DirectSumMat(m));
    od;
    d:=Length(g1[1]);
    if d<>Length(g2[1]) or d<>Length(p) or DeterminantMat(p)^2<>1 then
        return fail;
    else
        return List(g1,x->x^p)=g2;
    fi;
end;

StablyPermutationMCheckMat:= function(g,c1,c2,p)
    local h;
    h:=List(ConjugacyClassesSubgroups2(g),Representative);
    return StablyPermutationCheckHMat(g,h,c1,c2,p);
end;

StablyPermutationFCheckMat:= function(g,c1,c2,p)
    local tg,gg,d,th,mi,ms,o,h,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return StablyPermutationCheckHMat(g2,h2,c1,c2,p);
    fi;
end;

StablyPermutationCheckHGen:= function(g,hh,c1,c2)
    local gg,g1,g2,dx,m,i,j;
    gg:=List(hh,x->CosetRepresentation(g,x));
    Add(gg,GeneratorsOfGroup(g));
    g1:=[];
    g2:=[];
    for i in [1..Length(gg[1])] do
        m:=[];
        for j in [1..Length(gg)] do
            m:=Concatenation(m,List([1..c1[j]],x->gg[j][i]));
        od;
        Add(g1,DirectSumMat(m));
        m:=[];
        for j in [1..Length(gg)] do
            m:=Concatenation(m,List([1..c2[j]],x->gg[j][i]));
        od;
        Add(g2,DirectSumMat(m));
    od;
    return [g1,g2];
end;

StablyPermutationMCheckGen:= function(g,c1,c2)
    local h;
    h:=List(ConjugacyClassesSubgroups2(g),Representative);
    return StablyPermutationCheckHGen(g,h,c1,c2);
end;

StablyPermutationFCheckGen:= function(g,c1,c2)
    local tg,gg,d,th,mi,ms,o,h,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return StablyPermutationCheckHGen(g2,h2,c1,c2);
    fi;
end;

CheckCoflabbyResolutionBaseH:= function(g,hh,mi)
    local z0;
    z0:=List(hh,Z0lattice);
    return ForAll([1..Length(hh)],i->LatticeBasis(List(OrbitsDomain(hh[i],mi),Sum))=z0[i]);
end;

CheckCoflabbyResolutionBase:= function(g,mi)
    local hh,z0;
    hh:=List(ConjugacyClassesSubgroups2(g),Representative);
    return CheckCoflabbyResolutionBaseH(g,hh,mi);
end;

SearchCoflabbyResolutionBaseH:= function(g,hh,b)
    local z0,orbs,imgs,i,j,mi,mis;
    mis:=[];
    z0:=List(hh,Z0lattice);
    orbs:=Set(Union(z0),x->Orbit(g,x));
    imgs:=List(orbs,x->List(hh,y->LatticeBasis(List(OrbitsDomain(y,x),Sum))));
    if b=0 then
        for i in [1..Length(orbs)] do
            for j in Combinations([1..Length(orbs)],i) do
                if ForAll([1..Length(hh)],x->LatticeBasis(
                  Union(List(j,y->imgs[y][x])))=z0[x]) then
                    mi:=Union(List(j,x->orbs[x]));
                    if mis=[] or Length(mi)<Length(mis) then
                        mis:=mi;
                    fi;
                fi;
            od;
            if mis<>[] then
                return mis;
            fi;
       od;
    else
        for i in [1..b] do
            for j in Combinations([1..Length(orbs)],i) do
                if ForAll([1..Length(hh)],x->LatticeBasis(
                  Union(List(j,y->imgs[y][x])))=z0[x]) then
                    mi:=Union(List(j,x->orbs[x]));
                    Add(mis,mi);
                fi;
            od;
        od;
        return Set(mis);
    fi;
end;

SearchCoflabbyResolutionBase:= function(g,b)
    local hh;
    hh:=List(ConjugacyClassesSubgroups2(g),Representative);
    return SearchCoflabbyResolutionBaseH(g,hh,b);
end;

FlabbyResolutionFromBase:= function(g,mi)
    local tg,gg,d,th,ms,o,r,gg1,gg2,v1,mg,img;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    r:=Length(mi);
    o:=IdentityMat(r);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return rec(injection:=TransposedMat(mi),
                   surjection:=NullMat(r,0),
                   actionP:=TransposedMatrixGroup(Group(gg1,o))
        );
    else
        ms:=NullspaceIntMat(mi);
        v1:=NullspaceIntMat(TransposedMat(ms));
        mg:=Concatenation(v1,ms);
        img:=mg^-1;
        gg2:=List(gg1,x->mg*x*img);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        return rec(injection:=TransposedMat(mi),
                   surjection:=TransposedMat(ms),
                   actionP:=TransposedMatrixGroup(Group(gg1)),
                   actionF:=TransposedMatrixGroup(Group(gg2))
        );
    fi;
end;

PossibilityOfStablyPermutationFFromBase:= function(g,mi)
    local tg,gg,d,th,ms,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return [0,1];
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return PossibilityOfStablyPermutationH(g2,h2);
    fi;
end;

StablyPermutationFCheckFromBase:= function(g,mi,c1,c2)
    local tg,gg,d,th,ms,o,h,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return StablyPermutationCheckH(g2,h2,c1,c2);
    fi;
end;

StablyPermutationFCheckPFromBase:= function(g,mi,c1,c2)
    local tg,gg,d,th,ms,o,h,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return StablyPermutationCheckHP(g2,h2,c1,c2);
    fi;
end;

StablyPermutationFCheckMatFromBase:= function(g,mi,c1,c2,p)
    local tg,gg,d,th,ms,o,h,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return StablyPermutationCheckHMat(g2,h2,c1,c2,p);
    fi;
end;

StablyPermutationFCheckGenFromBase:= function(g,mi,c1,c2)
    local tg,gg,d,th,ms,o,h,r,gg1,gg2,g2,iso,ker,tg2,th2,h2,m1,m2,v;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=List(ConjugacyClassesSubgroups2(tg),Representative);
    r:=Length(mi);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return true;
    else
        ms:=NullspaceIntMat(mi);
        v:=NullspaceIntMat(TransposedMat(ms));
        m1:=Concatenation(v,ms);
        m2:=m1^-1;
        gg2:=List(gg1,x->m1*x*m2);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        tg2:=Group(gg2);
        iso:=GroupHomomorphismByImages(tg,tg2,gg,gg2);
        ker:=Kernel(iso);
        th2:=List(Filtered(th,x->IsSubgroup(x,ker)),x->Image(iso,x));
        g2:=TransposedMatrixGroup(tg2);
        h2:=List(th2,TransposedMatrixGroup);
        return StablyPermutationCheckHGen(g2,h2,c1,c2);
    fi;
end;

Norm1TorusJ :=function(d,n)
    local I,M1,M2,M,f,Sn,T;
    I:=IdentityMat(d-1);
    Sn:=SymmetricGroup(d);
    T:=TransitiveGroup(d,n);
    M1:=Concatenation(List([2..d-1],x->I[x]),[-List([1..d-1],One)]);
    if d=2 then
        M:=[M1];
    else
        M2:=Concatenation([I[2],I[1]],List([3..d-1],x->I[x]));
        M:=[M1,M2];
    fi;
    f:=GroupHomomorphismByImages(Sn,Group(M),GeneratorsOfGroup(Sn),M);
    return Image(f,T);
end;

DirectSumMatrixGroup:= function(l)
    local gg,gg1;
    gg:=List(l,GeneratorsOfGroup);
    if Length(Set(gg,Length))>1 then
        return fail;
    else
        gg1:=List([1..Length(gg[1])],x->DirectSumMat(List(gg,y->y[x])));
    fi;
    return Group(gg1,DirectSumMat(List(l,Identity)));
end;

DirectProductMatrixGroup:= function(l)
    local gg,gg1,o,o1,i,j,gx;
    gg:=List(l,GeneratorsOfGroup);
    gg1:=[];
    for i in [1..Length(l)] do
        o:=List(l,Identity);
        for j in gg[i] do
            o[i]:=j;
            Add(gg1,DirectSumMat(o));
        od;
    od;
    return Group(gg1,DirectSumMat(List(l,Identity)));
end;

FlabbyResolutionFromPerm:= function(g)
    local tg,gg,d,th,mi,ms,o,r,gg1,gg2,v1,mg,img;
    tg:=TransposedMatrixGroup(g);
    gg:=GeneratorsOfGroup(tg);
    d:=Length(Identity(g));
    th:=ConjugacyClassesSubgroupsFromPerm(tg);
    mi:=FindCoflabbyResolutionBase(tg,th);
    r:=Length(mi);
    o:=IdentityMat(r);
    gg1:=List(gg,x->PermutationMat(Permutation(x,mi),r));
    if r=d then
        return rec(injection:=TransposedMat(mi),
                   surjection:=NullMat(r,0),
                   actionP:=TransposedMatrixGroup(Group(gg1,o))
        );
    else
        ms:=NullspaceIntMat(mi);
        v1:=NullspaceIntMat(TransposedMat(ms));
        mg:=Concatenation(v1,ms);
        img:=mg^-1;
        gg2:=List(gg1,x->mg*x*img);
        gg2:=List(gg2,x->x{[d+1..r]}{[d+1..r]});
        return rec(injection:=TransposedMat(mi),
                   surjection:=TransposedMat(ms),
                   actionP:=TransposedMatrixGroup(Group(gg1)),
                   actionF:=TransposedMatrixGroup(Group(gg2))
        );
    fi;
end;

