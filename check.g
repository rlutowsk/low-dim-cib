CoboundaryBasisEquationByGens := function(gens)
    local equation, id, snf;

    id := One(gens[1]);
    equation := TransposedMat(
        Concatenation(
            List(
                gens,
                x->TransposedMat(x-id)
            )
        )
    );
    snf  := SmithNormalFormIntegerMatTransforms(equation);
    if snf.rank = 0 then
        return equation{[1]};
    fi;
    return (snf.coltrans^-1){[1..snf.rank]};
end;

CaratCatQName := function(str)
    local i,f,l; 
    l:=Length(str); 
    f:=false; 
    for i in [1..l] do
        if str[i]='.' then
            if f then break;
            else
                f:=true;
            fi;
        fi;
    od; 
    if f then return str{[1..i-1]}; fi; 
    return fail; 
end;

CaratCatZName := function(str)
    local i,f,l; 
    l:=Length(str); 
    f:=0; 
    for i in [1..l] do
        if str[i]='.' then
            f := f+1;
            if f = 4 then break; fi;
        fi;
    od; 
    if f=4 then return str{[1..i-1]}; fi; 
    return fail; 
end;

RedisCommandFormatted := function(arg) return RedisCommand(CallFuncList(StringFormatted,arg)); end;

CaratCatCoboundaryBasisEquationByZName := function( name )
    local zcl;
    if not IsString(name) then return fail; fi;
    zcl := RedisCommandFormatted( "HGET {} generators", name );
    if zcl = "" then return fail; fi;
    zcl := EvalString(zcl);
    return CoboundaryBasisEquationByGens( zcl );
end;

CheckCIBCohomology := function(aname, eqn...)
    local vs, cib, zcl, eq;

    vs  := RedisCommandFormatted("HGET {} vs", aname);
    if vs = "" then return fail; fi;
    cib := RedisCommandFormatted("HGET {} cib:f", aname);
    if cib = "" then return fail; fi;
    vs  := Concatenation( EvalString( vs ) );
    cib := List( CIBUnflatString( cib ), Concatenation );
    if IsBound(eqn[1]) and IsMatrix( eqn[1] ) then
        eq := eqn[1];
    else
        eq  := CaratCatCoboundaryBasisEquationByZName( CaratCatZName(aname) );
    fi;
    return List( cib, c->SolveInhomEquationsModZ( eq, c-vs, true )[1] <> [] );
end;

CheckCoboundary := function( zname )
    local z, s, cbz, cbs, mzs, msz;
    z := CaratCatFetchZClass( zname );
    s := CaratCatFetchAffineCrystGroup( Concatenation(zname,".0") );

    cbz := List( CIB.CoboundaryBasisIntOnRightByGens( z.generators ), Concatenation );
    cbs := List( CoboundaryBasisInt( s ), Concatenation );

    msz := ComplementIntMat( cbs, cbz );
    mzs := ComplementIntMat( cbz, cbs );

    return ForAll(msz.moduli, x->x=1) and ForAll(mzs.moduli, x->x=1);
end;

CaratCatRecalculate := function( aname )
    local s,r,c,d,l,f, use_data;

    if ValueOption("context") = true then
        s := CaratCatFetchAffineCrystGroup( aname );
        SetCofiniteIntegralBraceVectorSystems( s, CIB.CofiniteIntegralBraceVectorSystemsByContext(s) );
    else
        s := CaratCatCofiniteIntegralBraceVectorSystems( aname, "cib:f");
    fi;
    r := CofiniteIntegralBracesRepsGenerators( s );
    d := DegreeOfMatrixGroup( s );
    if r<>[] then
        l := Length( r[1] );
    else
        l := 0;
    fi;
    c:=CIBFlatString( List(r, x->x{[1..l]}[d]{[1..d-1]}) );

    f := RedisCommand("HGET {} cib:f", aname);
    if f="[]" then f:="[0,0,0,1]"; fi;
    return f=c;
    #if f <> c then
    #    Error("Values must match");
    #fi;
end;

CaratCatCheckSavedData := function( zname, prefix )
    local stream, data, anames, aff;

    Info( InfoCaratCat, 2, zname, ": start processing");
    
    stream := InputTextFile( Concatenation(prefix,"/",zname) );

    data := JsonStreamToGap( stream );

    CloseStream( stream );

    SortBy(data.aff, x->x.no);

    Info( InfoCaratCat, 2, zname, ": read json data");

    anames := SSortedList(RedisCommand("SMEMBERS {}:anames", zname), Int);
    if anames <> List(data.aff, x->x.no) then
        return false;
    fi;
    Info( InfoCaratCat, 2, zname, ": ", Size(anames), " affine classes match");
    for aff in data.aff do 
        if aff <> CaratCatAffRec2( zname, aff.no, "vs", "cib:f" ) then
            Info( InfoCaratCat, 2, zname, ": data mismatch");
            return false;
        fi;
    od;
    Info( InfoCaratCat, 1, zname, ": data match");
    return true;
end;
