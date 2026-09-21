CaratCatSaveCIBCocycles := function(sgrp, keys...)
    local cibs, dim, str, key;
    if not HasName(sgrp) then
        Error("The group must have a name. Use CaratCatFetchAffineCrystGroup function to get a group with the CARAT name");
    fi;

    if Length(keys) > 0 then
        key := keys[1];
    else
        key := "cib";
    fi;

    cibs := CofiniteIntegralBracesRepsGenerators( sgrp );
    Info(InfoCaratCat, 1, Name(sgrp), ": computed ", Size(cibs), " CIB representatives");
    if Length(cibs) = 0 then
        RedisCommand(StringFormatted("HSET {} {} []", Name(sgrp), key));
        return 0;
    fi;
    dim := DegreeOfMatrixGroup(sgrp)-1;
    if Size(PointGroup(sgrp))>1 and List(cibs[1], g->g{[1..dim]}{[1..dim]}) <> GeneratorsOfGroup( PointGroup(sgrp) ) then
        Error("The point group of the group must be the same as the point group of the Z-class. Use CaratCatFetchAffineCrystGroup function to get the group");
    fi;
    str := CIBFlatString( List(cibs, cib->List(cib, g->g[dim+1]{[1..dim]})) );
    RemoveCharacters(str, " ");
    RedisCommand(StringFormatted("HSET {} {} {}", Name(sgrp), key, str));
    return Size( cibs );
end;

CIBPrepareQueue := function(dim, queue, append...)
    local res, res1, x, y, i, n;
    if Length(append) > 0 and not append[1] then
        RedisCommand(StringFormatted("DEL {}", queue));
    fi;
    n := RedisCommand(StringFormatted("LLEN {}", queue));
    Info(InfoCaratCat, 1, "Redis list ", queue, " has ", n, " elements");
    res := RedisCommand(StringFormatted("SMEMBERS zclasses:dim:{}", dim));
    i := 0;
    for x in res do
        res1 := RedisCommand(StringFormatted("SMEMBERS {}:anames", x));
        for y in res1 do
             RedisCommand(StringFormatted("RPUSH {} {}.{}", queue, x, y));
             i := i+1;
        od;
    od;
    Info(InfoCaratCat, 1, "Invoked ", i, " times RPUSH to ", queue);
    return RedisCommand(StringFormatted("LLEN {}", queue)) - n - i;
end;

CIBFetchName := function(queue, set)
    local res;
    res := RedisCommand(StringFormatted("LPOP {}", queue));
    if res = "" then
        return fail;
    fi;
    RedisCommand(StringFormatted("SADD {} {}", set, res));
    return res;
end;

CIBRemoveName := function(name, set)
    RedisCommand(StringFormatted("SREM {} {}", set, name));
end;

CIBRunJob1 := function(queue, set, maxexp...)
    local sgrp, name, max, len, exp, exps;

    if Length(maxexp) > 0 then
        max := maxexp[1];
    else
        max := 0;
    fi;
    while true do
        name := CIBFetchName(queue, set);
        if name = fail then
            return;
        fi;
        sgrp := CaratCatFetchAffineCrystGroup(name);
        if IsSolvableGroup( sgrp ) then
            exps := BraceMaxExponents( sgrp );
        else
            exps := [];
        fi;
        SetBraceMaxExponents( sgrp, Immutable(exps) );
        if max > 0 and Length( BraceMaxExponents(sgrp) ) > 0 and Maximum( BraceMaxExponents(sgrp) ) > max then
            Info(InfoCaratCat, 1, name, ": reached threshold ", max, "; skipping ...");
            RedisCommand("MULTI");
            RedisCommand(StringFormatted("LLEN {}.skipped", queue));
            RedisCommand(StringFormatted("RPUSH {}.skipped {}", queue, name));
            len := RedisCommand("EXEC");
            if  len[2] <> len[1]+1 then
                Error("Failed to push back the <name> to the queue: ", name);
            fi;
            RedisCommand(StringFormatted("SREM {} {}", set, name));
            continue;
        fi;
        CaratCatSaveCIBCocycles(sgrp);
        RedisCommand(StringFormatted("SREM {} {}", set, name));
    od;
end;

# CIBBySubgroupPointGroup := function( sgrp, spg, cnt... )
#     local cb, ind, ssg, cib, vs, c, sol, cand, res, rk, m;

#     ind := List( SSortedList(spg), g -> Position( SSortedList( PointGroup(sgrp) ), g ) );
#     ssg := PreImage( PointHomomorphism(sgrp), spg );
#     if InfoLevel(InfoCIB) > 2 and Size(cnt)>0 and IsInt(cnt[1]) then
#         SetName(ssg, StringFormatted("{}:{}", Name(sgrp), cnt[1] ));
#         Info( InfoCIB, 3, Name(ssg), ": starting calculations");
#     fi;
#     if PointGroup(ssg) <> spg then
#         Error("The subgroup is not a subgroup of the point group of the group. This should not happen.");
#     fi;
#     vs  := VectorSystem( sgrp );
#     cib := List( CofiniteIntegralBraceVectorSystemsByContext( ssg ), x->Concatenation(x-vs{ind}));
#     # vs  := Concatenation( VectorSystem( sgrp ){ind} );
#     cb  := List( CoboundaryBasisInt( sgrp ), x -> Concatenation( x{ind} ) );
#     sol := [];
#     m   := Size( PointGroup(sgrp) );
#     for c in cib do
#         Append( sol, m * SolveInhomEquationsModZ(cb, c, true)[1]);
#     od;
#     # sol := Concatenation( List( cib, c -> SolveInhomEquationsModZ(cb, c, true)[1] ) );
#     res := [];
#     cb  := CoboundaryBasisInt( sgrp );
#     rk  := Size( cb );
#     vs  := vs * m;
#     for c in sol do
#         cand := GetZmodnZGroup( vs + Sum([1..rk], i->c[i]*cb[i]), m );
#         if cand<>fail then
#             Add( res, cand/m );
#         fi;
#     od;
#     return res;
# end;

# CIBByMaximalNormalSubgroupsPointGroup := function(sgrp, pp, checked)
#     local mns, res, spg;

#     Info( InfoCIB, 3, "Entering CIBByMaximalNormalSubgroupsPointGroup for size ", Size(PointGroup(sgrp)) );

#     mns := MaximalNormalSubgroups( PointGroup(sgrp) );
#     res := [];
#     if Size(mns[1]) > pp then
#         for spg in mns do
#             Append(res, CIBByMaximalNormalSubgroupsPointGroup( PreImage(PointHomomorphism(sgrp), spg), pp, checked ));
#         od;
#     else
#         for spg in mns do
#             if spg in checked then; continue; fi;
#             Add( checked, spg );
#             Append( res, CIBBySubgroupPointGroup(sgrp, spg) );
#         od;
#     fi;
#     return SSortedList( res );
# end;

# CIBBySylowSubgroupsPointGroup := function( sgrp )
#     local syl, res, p, pd, d, sg, iso, img;

#     if HasCofiniteIntegralBraceVectorSystems( sgrp ) then
#         return CofiniteIntegralBraceVectorSystems( sgrp );
#     fi;

#     if not IsSolvableGroup( sgrp ) then
#         return [];
#     fi;

#     pd := PrimeDivisors( Size(PointGroup(sgrp)) );
#     if pd = [] then
#         return ListWithIdenticalEntries( DegreeOfMatrixGroup(sgrp)-1, 0 );
#     fi;

#     res := [];
#     d   := 0;
#     if Size(pd) = 1 then
#         # handle p-groups
#         p := pd[1];
#         Info( InfoCIB, 1, Name(sgrp), ": calculating CIBs for ", p, "-group of order ", Size(PointGroup(sgrp)) );
#         if IsOddInt( p ) or Size(PointGroup(sgrp)) < 10 then
#             return CofiniteIntegralBraceVectorSystemsByContext( sgrp );
#         fi;
#         iso := IsomorphismPcGroup( PointGroup(sgrp) );
#         img := Image( iso );
#         for sg in Concatenation(
#             List(
#                 Filtered( 
#                     ConjugacyClassesSubgroups( img ), 
#                     x->Size(Representative(x))=8 
#                 ), 
#                 List)
#             ) do
#             d := d+1;
#             Append( res, CIBBySubgroupPointGroup( sgrp, PreImage(iso, sg), d ) );
#             res := Unique( res );
#         od;
#         return res;
#     fi;
#     # take first odd prime and calculate with this one
#     p := pd[2];
#     Info( InfoCIB, 1, Name(sgrp), ": calculating CIBs using ", p, "-Sylow subgroups of order ", Size(SylowSubgroup(PointGroup(sgrp), p)) );

#     for syl in SylowSubgroup( PointGroup( sgrp ), p )^PointGroup( sgrp ) do
#         d := d+1;
#         Append( res, CIBBySubgroupPointGroup( sgrp, syl, d ) );
#         res := Unique( res );
#     od;

#     return res;
# end;

CIBRunJob := function(queue, set)
    local sgrp, name, max, len, exp, cib;

    while true do
        name := CIBFetchName(queue, set);
        if name = fail then
            QuitGap(0);
        fi;
        # if RedisCommand(StringFormatted("HGET {} cib2", name)) <> "" then
        #     Info( InfoCaratCat, 1, name, ": already computed");
        #     continue;
        # fi;
        sgrp := CaratCatFetchAffineCrystGroup(name); # : subgroups:="cib1" );
        # if sgrp[2] = false then
        #     Info( InfoCaratCat, 1, name, ": found 0 CIB representatives by data");
        #     SetCofiniteIntegralBracesRepsGenerators( sgrp[1], [] );
        # fi;
        # sgrp := sgrp[1];
        # CofiniteIntegralBracesRepsGenerators( sgrp );
        cib := CIB.CofiniteIntegralBraceVectorSystemsByContext( sgrp );
        SetCofiniteIntegralBraceVectorSystems( sgrp, cib );
        # sgrp := CaratCatCofiniteIntegralBraceVectorSystems( name, "cib1" );
        CofiniteIntegralBracesRepsGenerators( sgrp );
        CaratCatSaveCIBCocycles(sgrp, "cib" );
        RedisCommand(StringFormatted("SREM {} {}", set, name));
        Unbind(sgrp);
    od;
end;

CIBRunJob1 := function(queue, set)
    local sgrp, name, max, len, exp, cib;

    while true do
        name := CIBFetchName(queue, set);
        if name = fail then
            QuitGap(0);
        fi;
        # if RedisCommand(StringFormatted("HGET {} cib2", name)) <> "" then
        #     Info( InfoCaratCat, 1, name, ": already computed");
        #     continue;
        # fi;
        # sgrp := CaratCatFetchAffineCrystGroup(name); # : subgroups:="cib1" );
        # if sgrp[2] = false then
        #     Info( InfoCaratCat, 1, name, ": found 0 CIB representatives by data");
        #     SetCofiniteIntegralBracesRepsGenerators( sgrp[1], [] );
        # fi;
        # sgrp := sgrp[1];
        # CofiniteIntegralBracesRepsGenerators( sgrp );
        # cib := CIB.CofiniteIntegralBraceVectorSystemsByContext( sgrp );
        # SetCofiniteIntegralBraceVectorSystems( sgrp, cib );
        sgrp := CaratCatCofiniteIntegralBraceVectorSystems( name, "cib" );
        CofiniteIntegralBracesRepsGenerators( sgrp );
        CaratCatSaveCIBCocycles(sgrp, "cib1" );
        RedisCommand(StringFormatted("SREM {} {}", set, name));
        Unbind(sgrp);
    od;
end;
