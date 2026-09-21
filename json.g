LoadPackage("json");

CaratCatAffRec := function(zname, aname, vskey, cibkey)
    local list, i, len, res;

    list := RedisCommand("HGETALL {}.{}", zname, aname);
    len  := Length( list );

    res := rec(no := Int(aname) );
    for i in [1,3..len-1] do
        if list[i] = vskey then
            res.vs := CIBFlatString( [ EvalString(list[i+1]) ] );
        elif list[i] = cibkey then
            res.cib := list[i+1];
        fi;
    od;
    return res;
end;

CaratCatZToJson := function( zname, vskey, cibkey )
    local res, a;
    res := CaratCatFetchZClass( zname );

    res.generators := CIBFlatString( res.generators );
    res.normalizer := CIBFlatString( res.normalizer );

    res.aff := [];
    for a in SSortedList(RedisCommand("SMEMBERS {}:anames", zname)) do
        Add( res.aff, CaratCatAffRec( zname, a, vskey, cibkey ) );
    od;

    return GapToJsonString( res );
end;

CaratCatZToJsonStream := function( zname, vskey, cibkey, stream )
    local res, a;
    res := CaratCatFetchZClass( zname );

    res.generators := CIBFlatString( res.generators );
    res.normalizer := CIBFlatString( res.normalizer );

    res.aff := [];
    for a in SSortedList(RedisCommand("SMEMBERS {}:anames", zname)) do
        Add( res.aff, CaratCatAffRec( zname, a, vskey, cibkey ) );
    od;

    GapToJsonStream( stream, res );
end;

CaratCatZToJsonGZ := function( zname, vskey, cibkey, prefix )
    local stream;

    stream := OutputGzipFile( Concatenation(prefix,"/",zname,".gz"), false );

    CaratCatZToJsonStream( zname, vskey, cibkey, stream );

    CloseStream( stream );
end;

CaratCatJsonGZToZ := function( zname, prefix )
    local stream, res, x;

    stream := InputTextFile( Concatenation(prefix,"/",zname) );

    res := JsonStreamToGap( stream );

    CloseStream( stream );

    for x in res.aff do
        if IsBound(x.cib) then
            x.cib := CIBUnflatString(x.cib);
        else
            x.cib := "";
        fi;
        x.vs  := CIBUnflatString(x.vs)[1];
    od;
    res.generators := CIBUnflatString(res.generators);
    res.normalizer := CIBUnflatString(res.normalizer);

    return res;
end;

CaratCatAffRec1 := function(zname, aname, vskey, cibkey)
    local list, i, len, res;

    list := RedisCommand("HGETALL {}.{}", zname, aname);
    len  := Length( list );

    res := rec(no := Int(aname) );
    for i in [1,3..len-1] do
        if list[i] = vskey then
            res.vs := EvalString(list[i+1]);
        elif list[i] = cibkey then
            res.cib := CIBUnflatString(list[i+1]);
        fi;
    od;
    if not IsBound(res.cib) then
        res.cib := "";
    fi;
    return res;
end;

CaratCatAffRec2 := function(zname, aname, vskey, cibkey)
    local list, i, len, res;

    list := RedisCommand("HGETALL {}.{}", zname, aname);
    len  := Length( list );

    res := rec(no := Int(aname) );
    for i in [1,3..len-1] do
        if list[i] = vskey then
            res.vs := CIBFlatString([EvalString(list[i+1])]);
        elif list[i] = cibkey then
            res.cib := list[i+1];
        fi;
    od;
    if not IsBound(res.cib) then
        res.cib := "";
    fi;
    return res;
end;


CaratCatFetchZClassFull := function( zname, keys... )
    local res, a, vskey, cibkey;
    vskey := "vs";
    cibkey:= "cib:f";
    if IsList(keys) then
        if IsBound(keys[1]) then
            vskey := keys[1];
        fi;
        if IsBound(keys[2]) then
            cibkey := keys[2];
        fi;
    fi;
    res := CaratCatFetchZClass( zname );
    res.aff := [];
    for a in SSortedList(RedisCommand("SMEMBERS {}:anames", zname)) do
        Add( res.aff, CaratCatAffRec1( zname, a, vskey, cibkey ) );
    od;
    return res;
end;
