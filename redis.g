CIBFlatString := function(strin)
    local res, flat, denom, str, list;

    if strin = [] then
        return "[0,0,0,1]";
    fi;
    if IsString(strin) then
        list   := EvalString( strin );
    else
        list   := strin;
    fi;

    res    := [];
    res[1] := Length( list );
    res[2] := Length( list[1] );
    res[3] := Length( list[1,1] );
    flat   := Flat( list );
    denom  := Lcm( List( flat, DenominatorRat ) );
    res[4] := denom;

    Append( res, denom * flat );

    str := String( res );
    RemoveCharacters( str, " " );

    return str;
end;

CIBUnflatString := function( strin )
    local list, num_mats, rows, cols, denom, flat_ints, 
          flat_rats, res, i, j, k, idx;
    
    if strin = "[]" then
        return [];
    fi;

    list := EvalString( strin );

    num_mats := list[1];
    rows     := list[2];
    cols     := list[3];
    denom    := list[4];

    flat_rats := list{[5..Length(list)]} / denom;

    res := [];
    idx := 1;
    for i in [1 .. num_mats] do
        res[i] := [];
        for j in [1 .. rows] do
            res[i][j] := [];
            for k in [1 .. cols] do
                res[i][j][k] := flat_rats[idx];
                idx := idx + 1;
            od;
        od;
    od;
    
    return res;
end;

if RedisConnected() then
    RedisFree();
fi;

RedisConnect("localhost", 6379);

CaratCatSaveCocycles := function()
    local cmd1, cmd2, x, i, siz, str, len;
    if not IsRecord(CDB) or not IsBound(CDB.root_dir) then
        Error("CDB record not initialized. Run CaratCatInit() first.");
    fi;

    # We work only with the right action
    if not IsBound(CDB.CocyclesRight) then
        Error("CDB.CocyclesRight not initialized. Compute or read the data first.");
    fi;
    i := 0;
    siz := Size(CDB.CocyclesRight);
    len := Length(String(siz));
    for x in CDB.CocyclesRight do
        str := String(x.mat);
        RemoveCharacters(str, " ");
        # Add data to the has set
        cmd1 := StringFormatted("HSET {}.{} vs {} orbit {}", x.zname, x.name, str, x.orbit);
        RedisCommand(cmd1);
        # Add info about the group to the set
        cmd2 := StringFormatted("SADD {}:anames {}", x.zname, x.name);
        RedisCommand(cmd2);
        if i mod 1000000 = 0 then
            Info(InfoCaratCat, 2, "[", PrintString(i, len), "/", siz, "]: ", cmd1);
            Info(InfoCaratCat, 2, "[", PrintString(i, len), "/", siz, "]: ", cmd2);
        fi;
        i := i + 1;
    od;
end;

CaratCatSaveZClasses := function()
    local cmd, x, i, siz, str, len;
    if not IsRecord(CDB) or not IsBound(CDB.root_dir) then
        Error("CDB record not initialized. Run CaratCatInit() first.");
    fi;

    # We work only with the right action
    if not IsBound(CDB.ZClassesRight) then
        Error("CDB.ZClassesRight not initialized. Compute or read the data first.");
    fi;
    i := 0;
    siz := Size(CDB.ZClassesRight);
    len := Length(String(siz));
    for x in CDB.ZClassesRight do
        if not IsBound(x.size) then
            x.size := Size( Group( x.generators ) );
        fi;
        str := String(x.generators);
        RemoveCharacters(str, " ");
        cmd := StringFormatted("HSET {} generators {} size {} dim {}", x.name, str, x.size, x.dim);
        RedisCommand(cmd);

        str := String(x.normalizer);
        RemoveCharacters(str, " ");
        cmd := StringFormatted("HSET {} normalizer {}", x.name, str);
        RedisCommand(cmd);

        str := String(x.formspace);
        RemoveCharacters(str, " ");
        cmd := StringFormatted("HSET {} formspace {}", x.name, str);
        RedisCommand(cmd);

        cmd := StringFormatted("SADD zclasses:dim:{} {}", x.dim, x.name);
        RedisCommand(cmd);

        if i mod 10000 = 0 then
            Info(InfoCaratCat, 2, "[", PrintString(i, len), "/", siz, "]: Processed ", x.name);
        fi;
        i := i + 1;
    od;
end;

CaratCatSaveQClasses := function(dim)
    local zclasses, z, s;

    zclasses := RedisCommand(StringFormatted("SMEMBERS zclasses:dim:{}", dim));
    for z in zclasses do
        s := SplitString(z, ".");
        RedisCommand(StringFormatted("SADD {}.{} {}.{}", s[1], s[2], s[3], s[4]));
    od;
end;

CaratCatSaveTorsionfree := function()
    local x;
    if not IsRecord(CDB) or not IsBound(CDB.root_dir) then
        Error("CDB record not initialized. Run CaratCatInit() first.");
    fi;

    # We work only with the right action
    if not IsBound(CDB.Torsionfree) then
        Error("CDB.Torsionfree not initialized. Compute or read the data first.");
    fi;
    for x in CDB.Torsionfree do
        RedisCommand(StringFormatted("SADD torsionfree {}", x));
    od;
end;

CaratCatFetchTorsionfree := function()
    return RedisCommand(StringFormatted("SMEMBERS torsionfree"));
end;

CaratCatFetchZClass := function(name)
    local res, rval;

    res := RedisCommand(StringFormatted("HGET {} generators", name));
    if res = "" then
        return fail;
    fi;
    rval := rec( generators := EvalString(res) );
    res := RedisCommand(StringFormatted("HGET {} normalizer", name));
    rval.normalizer := EvalString(res);
    res := RedisCommand(StringFormatted("HGET {} size", name));
    rval.size := Int(res);
    if rval.size = 0 then
        Error("Size of the group cannot be zero. Data might be corrupted for ", name);
    fi;
    res := RedisCommand(StringFormatted("HGET {} dim", name));
    rval.dim := Int(res);
    if rval.dim = 0 then
        Error("Dimension of the group cannot be zero. Data might be corrupted for ", name);
    fi;
    rval.name := name;
    return rval;
end;

MaximalPSubgroups := function( grp )
    local pd;

    pd := PrimeDivisors( Size(grp) );
    if Size(pd)=1 then
        return MaximalNormalSubgroups( grp );
    fi;
    return List(pd, p->SylowSubgroup(grp,p));
end;

CaratCatZClassCatalog := function( grp )
    local res, names, n, grp1, zcl;

    res   := CaratQClassCatalog( grp, 0 );
    names := RedisCommand(StringFormatted("SMEMBERS {}", res.qclass));
    for n in names do
        zcl := CaratCatFetchZClass(Concatenation(res.qclass,".",n));
        grp1:= Group( zcl.generators );
        if RepresentativeAction( GL(DegreeOfMatrixGroup(grp), Integers), grp1, grp ) <> fail then
            return zcl.name;
        fi;
    od;
    return fail;
end;

CaratCatSaveMaxPSubgroups := function(name)
    local r, pg, grp, zname, str, list;

    if RedisCommand(StringFormatted("HEXISTS {} maxpsubs", name)) = 1 then
        return;
    fi;

    r := CaratCatFetchZClass( name );
    pg:= Group( r.generators );

    list := [];
    for grp in MaximalPSubgroups( pg ) do
        zname := CaratCatZClassCatalog( grp );
        # zname := CaratName( grp );
        Add( list, zname );
    od;
    str := String(list);
    RemoveCharacters(str, " ");
    RedisCommand(StringFormatted("HSET {} maxpsubs {}", name, str));
end;

CaratCatFetchCocycle := function(name)
    local res;

    res := RedisCommand(StringFormatted("HGET {} vs", name));
    if res = "" then
        return fail;
    fi;
    return EvalString(res);
end;

VectorSystemByGensOnRight := function( sgens )
    local orbit_s, orbit_p, d, p, i, gens, nd, g, dd, n, nn, exp;

    d  := Size( sgens[1] ) - 1;
    nd := [1..d];
    orbit_s := [ IdentityMat(d+1) ];
    orbit_p := [ IdentityMat(d)   ];

    gens := List( Filtered( sgens, x->x{nd}{nd}<>orbit_p[1] ), MutableCopyMat );
    exp  := Lcm( Concatenation( List( gens, s->List(s[d+1]{nd}, DenominatorRat) ) ) );
    for g in gens do
        g[d+1]{nd} := g[d+1]{nd} * exp;
    od;

    i := 1;
    while i <= Length( orbit_p ) do
        dd := orbit_s[i];
        for g in gens do
            n  := dd * g;
            nn := n{nd}{nd};
            if not nn in orbit_p then
                Add(orbit_p, nn);
                Add(orbit_s, n);
            fi;
        od;
        i := i+1;
    od;
    SortParallel( orbit_p, orbit_s );

    return ( ( List( orbit_s, x->x[d+1]{nd} ) ) mod exp ) / exp;
end;

CaratCatAdjustCIBData := function( sgrp, key )
    local split, zname, src, conj, P, P1, Plist, cib_str, cib, Sgens, ngens, c, d, nd, d1, ng, vs, N, Ngens, Npair, Nperm, nr, act, Adata, Psize;
    
    split := SplitString(Name(sgrp), ".");
    zname := Concatenation(split[1], ".", split[2], ".", split[3], ".", split[4]);
    src   := CaratCatFetchZClass( zname );

    P     := PointGroup( sgrp );
    Plist := SSortedList( P );
    Psize := Size( P );
    P1    := Group(src.generators);
    conj  := RepresentativeAction( GL(DegreeOfMatrixGroup(P), Integers), P1, P );

    if P1^conj <> P then
        Error("Mismatch in calculation of conjugation matrix from PointGroup(<sgrp>) to CARAT group");
    fi;
    Sgens := List( src.generators, x->MutableCopyMat( DirectSumMat(x^conj, [[1]]) ) );
    ngens := Size( Sgens );
    
    cib_str := RedisCommand(StringFormatted("HGET {} {}", Name(sgrp), key));
    
    cib   := CIBUnflatString(cib_str);

    d     := DegreeOfMatrixGroup( P );
    ng    := [1..ngens];
    nd    := [1..d];
    d1    := d+1;
    vs    := [];
    
    N     := Group( Concatenation(src.generators, src.normalizer) )^conj;

    Ngens := [ One(N) ];
    Append( Ngens, GeneratorsOfGroup( N ) );
    Npair := List( Ngens, n->[n^-1, n] );
    Nperm := List( Npair, n->PermList( List(Plist, e->PositionSet( Plist, n[1] * e * n[2] )) ) );
    nr    := [1..Psize];
    act := function( list, pos )
        local img;
        img := Permuted( list, Nperm[pos] );
        img := MutableCopyMat( img * Ngens[pos] ) mod Psize; 
        return img;
    end;
    for c in cib do
        if Size(c) <> ngens then
            Error("Length of cib vector does not match the size of the generating set for <sgrp>");
        fi;
        Sgens{ng}[d1]{nd} := c * conj;
        Append( vs, Orbit(N, VectorSystemByGensOnRight(Sgens) * Psize, Ngens{[2..Size(Ngens)]}, [2..Size(Ngens)], act)/Psize );
    od;

    return rec( grp := sgrp, conj := conj, cib := vs );
end;

CaratCatAppendSubgroupsInfoP  := function(sgrp, key)
    local P, Psubs, Ssubs, s, res, iso;

    P := PointGroup( sgrp );
    if not IsPGroup( P ) then
        Error("Point group of <sgrp> must be a p-group");
    fi;
    iso   := IsomorphismPcGroup( P );
    Psubs := List( MaximalNormalSubgroups( Image(iso) ), s->PreImage(iso,s) );
    Ssubs := List( Psubs, x->PreImage(PointHomomorphism(sgrp),x) );
    res   := [];
    for s in Ssubs do
        SetName( s, CaratName(s) );
        Add(res, CaratCatAdjustCIBData(s,key));
    od;
    return [ sgrp, res ];
end;

MaxPrimePowerInt := function( num )
    local ppi, lp, lv;
    
    ppi := PrimePowersInt( num );
    lp  := ppi{[1,3..Size(ppi)-1]};
    lv  := List([1,3..Size(ppi)-1], i->ppi[i]^ppi[i+1]);
    
    return lp[Position(lv, Maximum(lv))];
end;

CaratCatAppendSubgroupsInfoPQ := function(sgrp, key)
    local P, Psyl, Ssubs, s, res, pd, p, name;

    P    := PointGroup( sgrp );
    pd   := PrimePowersInt( Size(P) );
    if Size(pd) < 2 then
        Error("This method is not for p-groups");
    fi;
    p    := MaxPrimePowerInt( Size(P) );
    Info( InfoCaratCat, 3, "Choosing prime ", p);
    Psyl := List( SylowSubgroup(P, p)^P );
    Ssubs:= List( Psyl, x->PreImage(PointHomomorphism(sgrp),x) );
    res  := [];
    name := CaratName( Ssubs[1] );
    for s in Ssubs do
        SetName( s, name );
        Add( res, CaratCatAdjustCIBData(s,key) );
    od;
    return [ sgrp, res ];
end;

CaratCatAppendSubgroupsInfo := function( sgrp, key )
    local P;
    P := PointGroup( sgrp );
    if not IsSolvable( P ) or Size( P )=1 then
        return [ sgrp, [] ];
    fi;
    if Size( PrimeDivisors( Size( P ) ) ) = 1 then
        return CaratCatAppendSubgroupsInfoP( sgrp, key );
    fi;
    return CaratCatAppendSubgroupsInfoPQ( sgrp, key );
end;

CaratCatBraceMaxExpAffineCrystGroup := function(name)
    local exps, split, qname, zname, z, g;
    
    split := SplitString(name, ".");
    qname := Concatenation(split[1], ".", split[2]);

    exps  := RedisCommand( StringFormatted("HGET exponents {}", qname) );
    if exps = "" then
        Info( InfoCIB, 1, qname, ": calculating max exponents" );
        zname := Concatenation(split[1], ".", split[2], ".", split[3], ".", split[4]);
        z := CaratCatFetchZClass(zname);
        g := Group(z.generators); IsGroup(g); IsFinite(g);
        exps := String( BraceMaxExponents( g ) );
        RemoveCharacters( exps, " " );
        RedisCommand( StringFormatted("HSET exponents {} {}", qname, exps) );
    fi;
end;

CaratCatFetchAffineCrystGroup := function(name)
    local z, c, zname, qname, split, sgens, i, mat, sgrp, P, N, exps, subgroups;

    split := SplitString(name, ".");
    zname := Concatenation(split[1], ".", split[2], ".", split[3], ".", split[4]);
    qname := Concatenation(split[1], ".", split[2]);
    z := CaratCatFetchZClass(zname);
    c := CaratCatFetchCocycle(name);
    sgens := [];
    for i in [1..Length(z.generators)] do
        mat := NullMat( z.dim+1, z.dim+1);
        mat{[1..z.dim]}{[1..z.dim]} := z.generators[i];
        mat[z.dim+1]{[1..z.dim]} := c[i];
        mat[z.dim+1, z.dim+1] := 1;
        Add(sgens, mat);
    od;
    sgrp := AffineCrystGroupOnRight(Concatenation(sgens, GeneratorsOfIntegralAffineCrystGroupOnRight(sgens)));

    SetIsStandardAffineCrystGroup(sgrp, true);

    P := Group(z.generators);
    N := Group(Concatenation(z.generators, z.normalizer));

    SetPointGroup(sgrp, P);

    P := PointGroup( sgrp );
    SetNormalizerInGLnZ(P, N);
    
    SetName(sgrp, name);

    if ValueOption("calculate") = true then
        exps := RedisCommand( StringFormatted("HGET exponents {}", qname) );
        if exps = "" then
            Info( InfoCIB, 1, qname, ": calculating max exponents" );
            exps := String( BraceMaxExponents( sgrp ) );
            RemoveCharacters( exps, " " );
            RedisCommand( StringFormatted("HSET exponents {} {}", qname, exps) );
        else
            exps := EvalString( exps );
            SetBraceMaxExponents( sgrp, exps );
            SetBraceMaxExponents( P, exps );
        fi;
    fi;

    subgroups := ValueOption("subgroups");
    if subgroups = "" then
        subgroups := false;
    fi;
    if not IsString(subgroups) then
        return sgrp;
    fi;
    return CaratCatAppendSubgroupsInfo( sgrp, subgroups );
end;

CaratCatCofiniteIntegralBraceVectorSystems := function( name, key )
    local data, sgrp, cibs, r, ssg, spg;

    Info( InfoCaratCat, 2, "Looking up database for subgroups of: ", name);
    
    data := CaratCatFetchAffineCrystGroup( name : subgroups:=key );
    sgrp := data[1];
    
    cibs := [];
    for r in data[2] do
        ssg := r.grp;
        Info(InfoCaratCat, 2, "Extending ", Name(ssg));
        spg := PointGroup( ssg );
        Append( cibs, CIB.ExtendedVectorSystems( sgrp, ssg, spg, r.cib ));
    od;
    
    cibs := Set( cibs );
    
    SetCofiniteIntegralBraceVectorSystems( sgrp, cibs );
    
    return sgrp;
end;

CaratCatCheckCofiniteIntegralBraceVectorSystems := function( name, key )
    local data, sgrp, cibs, r, ssg, spg, orig, cib1, d;

    Info( InfoCaratCat, 2, "Looking up database for subgroups of: ", name);
    
    data := CaratCatFetchAffineCrystGroup( name : subgroups:=key );
    sgrp := data[1];
    
    d := DegreeOfMatrixGroup(sgrp)-1;
    
    cibs := [];
    for r in data[2] do
        ssg := r.grp;
        Info(InfoCaratCat, 2, "Extending ", Name(ssg));
        spg := PointGroup( ssg );
        Append( cibs, CIB.ExtendedVectorSystems( sgrp, ssg, spg, r.cib ));
    od;
    
    cibs := Set( cibs );
    SetCofiniteIntegralBraceVectorSystems( sgrp, cibs );
    cibs := List( CofiniteIntegralBracesRepsGenerators( sgrp ), gens->List(gens, g->g[d+1]{[1..d]}) );

    orig := RedisCommand(StringFormatted("HGET {} {}", Name(sgrp), key));

    #sgrp := CaratCatFetchAffineCrystGroup(name);
    #cib1 := CIB.CofiniteIntegralBraceVectorSystemsByContext( sgrp );
    #SetCofiniteIntegralBraceVectorSystems( sgrp, cib1 );
    #cib1 := List( CofiniteIntegralBracesRepsGenerators( sgrp ), gens->List(gens, g->g[d+1]{[1..d]}) );
    
    return rec( sgrp:=sgrp, orig:=orig, calc:=CIBFlatString(cibs) );
end;
