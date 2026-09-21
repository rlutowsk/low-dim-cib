if not (IsBound(InfoCaratCat) and IsInfoClass(InfoCaratCat)) then
    InfoCaratCat := NewInfoClass("InfoCaratCat");
fi;

if not (IsBound(CDB) and IsRecord(CDB)) then
    BindGlobal("CDB", rec());
fi;

CaratReadMetadata := function(str)
    local comment, res, sub, pos, i;

    res := rec();
    comment := Position(str, '%');
    if comment = fail then
        return res;
    fi;
    sub := SplitString(str{[comment..Length(str)]}, ":, \n");
    pos := Position(sub, "orbit");
    if pos <> fail then
        res.orbit := Int(sub[pos+1]);
    fi;
    pos := Position(sub, "name");
    if pos <> fail then
        for i in [pos+1..Size(sub)] do
            if sub[i] = "" then
                continue;
            fi;
            res.name := Int(sub[i]);
            if IsInt(res.name) then
                break;
            else
                Error("Should be integer");
            fi;
        od;
    fi;
    return res;
end;

CaratReadMatricesMetadata := function ( input, n )
    local res, i, str, r, pos;
    res := [  ];
    for i in [ 1 .. n ] do
        str := ReadLine( input );
        r   := CaratReadMetadata( str );
        pos := Position(str, '%');
        if pos <> fail then
            str := str{[1..pos-1]};
        fi;
        r.mat := CaratReadMatrix( input, str );
        Add( res, r );
    od;
    return res;
end;

CaratReadMatrixFileMetadata := function ( filename )
    local input, str, pos, n, res, i;
    input := InputTextFile( filename );
    str := CaratReadLine( input );
    pos := Position( str, '#' );
    if pos <> fail then
        n := CaratNextNumber( str, pos + 1 );
        res := CaratReadMatricesMetadata( input, n );
    else
        res := CaratReadMatrix( input, str );
    fi;
    CloseStream( input );
    return res;
end;

CaratCatInit := function(dir)
    if not IsDirectoryPath( dir ) then
        Error("<dir> must point to a directory");
    fi;
    CDB.root_dir := dir;
    Info(InfoCaratCat, 1, "Set root directory to ", dir);
end;

CaratCatReadCocycleData := function()
    local root_dir, orig_dir, dim, work_dir, x, y, r, sizes, lengths, i, filename;

    if not IsRecord(CDB) or not IsBound(CDB.root_dir) then
        Error("CDB record not initialized. Run CaratCatInit() first.");
    fi;

    CDB.CocyclesLeft := [];

    sizes  := [2, 13, 73, 710, 6079, 85308];
    lengths:= [1,  2,  2,   3,    4,     5];

    for dim in [1..6] do
        Info(InfoCaratCat, 1, "Processing dimension ", dim, " ...");
        work_dir := Concatenation(CDB.root_dir, "/", String(dim),"/c");
        if not IsDirectoryPath(work_dir) then
            Error("Directory must exist: ", work_dir);
        fi;
        i := 0;
        for x in DirectoryContents(work_dir) do
            filename := Concatenation(work_dir,"/",x);
            if IsDirectoryPath(filename) then
                continue;
            fi;
            i := i+1;
            Info(InfoCaratCat, 2, "[", PrintString(i, lengths[dim]), "/", sizes[dim], "] Processing ", x, " ...");
            r := CaratReadMatrixFileMetadata(filename);
            for y in r do
                y.zname := x;
                y.dim   := dim;
                if not IsBound(y.name) then
                    y.name := 0;
                    y.orbit:= 1;
                fi;
            od;
            Append(CDB.CocyclesLeft, r);
        od;
    od;
end;

CaratCatReadZClassData := function()
    local dim, work_dir, x, y, r, sizes, lengths, i, filename;

    CDB.ZClassesLeft := [];

    sizes  := [2, 13, 73, 710, 6079, 85308];
    lengths:= [1,  2,  2,   3,    4,     5];

    for dim in [1..6] do
        Info(InfoCaratCat, 1, "Processing dimension ", dim, " ...");
        work_dir := Concatenation(CDB.root_dir, "/", String(dim),"/z");
        if not IsDirectoryPath(work_dir) then
            Error("Directory must exist: ", work_dir);
        fi;
        i := 0;
        for x in DirectoryContents(work_dir) do
            filename := Concatenation(work_dir,"/",x);
            if IsDirectoryPath(filename) then
                continue;
            fi;
            i := i+1;
            Info(InfoCaratCat, 2, "[", PrintString(i, lengths[dim]), "/", sizes[dim], "] Processing ", x, " ...");
            r := CaratReadBravaisFile(filename);
            r.name := x;
            r.dim  := dim;
            Add(CDB.ZClassesLeft, r);
        od;
    od;
end;

CaratCatNormalizeZClassData := function()
    local i, x, siz, N, len, skip, qcl;

    siz := Size(CDB.ZClassesLeft);
    len := Length( String(siz) );
    skip := [];
    for i in [1..siz] do
        x := CDB.ZClassesLeft[i];
        Info(InfoCaratCat, 2, "[", PrintString(i, len), "/", siz, "] Processing ", x.name, " ...");
        x.normalizer := Concatenation(x.generators, x.normalizer);

        qcl := SplitString(x.name, ".");
        qcl := Concatenation(qcl[1], ".", qcl[2]);
        if qcl in skip then
            continue;
        fi;

        N := Group( x.normalizer );
        if IsFinite(N) then
            x.normalizer := GeneratorsSmallest( N );
            Info(InfoCaratCat, 1, "[", PrintString(i, len), "/", siz, "] Normalized ", x.name);
        else
            Add(skip, qcl);
        fi;
    od;
end;

CaratCatWriteCocycleData := function()
    local i, str, outl, outr, x, siz, r, strl, strr, soutl, soutr;

    outl := OutputGzipFile( Concatenation( CDB.root_dir, "/c.left.grp.gz" ), false );
    outr := OutputGzipFile( Concatenation( CDB.root_dir, "/c.right.grp.gz" ), false );

    strl := "";
    strr := "";

    AppendTo(outl, "CDB.CocyclesLeft := [");
    AppendTo(outr, "CDB.CocyclesRight := [");
    i := 0;
    siz := Size(CDB.CocyclesLeft);

    for x in CDB.CocyclesLeft do
        if i mod 100000 = 0 then
            if i>0 then
                CloseStream(soutl);
                CloseStream(soutr);
                Info(InfoCaratCat, 2, "Writing batch no ", i/100000 - 1, " ...");
                AppendTo(outl, strl);
                AppendTo(outr, strr);
            fi;    
            soutl := OutputTextString(strl, false);
            soutr := OutputTextString(strr, false);
            Info(InfoCaratCat, 2, "Starting batch no ", i/100000, " ...");
        fi;
        str := String(x);
        RemoveCharacters(str, " ");
        AppendTo(soutl, "\n", str);

        x.mat := TransposedMat( x.mat );
        str := String(x);
        RemoveCharacters(str, " ");
        AppendTo(soutr, "\n", str);
        i := i+1;
        if i<siz then
            AppendTo(soutl,",");
            AppendTo(soutr,",");
        fi;
    od;
    if not IsClosedStream(soutl) then
        CloseStream(soutl);
        AppendTo(outl, strl);
    fi;
    if not IsClosedStream(soutr) then
        CloseStream(soutr);
        AppendTo(outr, strr);
    fi;
    AppendTo(outl, "\n];");
    AppendTo(outr, "\n];");
    CloseStream( outl );
    CloseStream( outr );
end;

CaratCatWriteZClassData := function()
    local i, str, outl, outr, x, siz, r, strl, strr, soutl, soutr, filename;

    if not ( IsRecord(CDB) and IsBound(CDB.root_dir) ) then
        Error("CDB record not initialized. Run CaratCatInit() first.");
    fi;

    filename := Concatenation( CDB.root_dir, "/z.left.grp.gz" );
    outl := OutputGzipFile( filename , false );
    Info( InfoCaratCat, 2, "Writing left action to ", filename, " ...");
    filename := Concatenation( CDB.root_dir, "/z.right.grp.gz" );
    outr := OutputGzipFile( filename , false );
    Info( InfoCaratCat, 2, "Writing right action to ", filename, " ...");

    strl := "";
    strr := "";

    AppendTo(outl, "CDB.ZClassesLeft := [");
    AppendTo(outr, "CDB.ZClassesRight := [");
    i := 0;
    siz := Size(CDB.ZClassesLeft);

    for x in CDB.ZClassesLeft do
        if i mod 10000 = 0 then
            if i>0 then
                CloseStream(soutl);
                CloseStream(soutr);
                Info(InfoCaratCat, 2, "Writing batch no ", i/10000 - 1, " ...");
                AppendTo(outl, strl);
                AppendTo(outr, strr);
            fi;    
            soutl := OutputTextString(strl, false);
            soutr := OutputTextString(strr, false);
            Info(InfoCaratCat, 2, "Starting batch no ", i/10000, " ...");
        fi;
        str := String(x);
        RemoveCharacters(str, " ");
        AppendTo(soutl, "\n", str);

        x.generators := List(x.generators, TransposedMat);
        x.normalizer := List(x.normalizer, TransposedMat);
        x.formspace  := List(x.formspace, TransposedMat);
        str := String(x);
        RemoveCharacters(str, " ");
        AppendTo(soutr, "\n", str);
        i := i+1;
        if i<siz then
            AppendTo(soutl,",");
            AppendTo(soutr,",");
        fi;
    od;
    if not IsClosedStream(soutl) then
        CloseStream(soutl);
        AppendTo(outl, strl);
    fi;
    if not IsClosedStream(soutr) then
        CloseStream(soutr);
        AppendTo(outr, strr);
    fi;
    AppendTo(outl, "\n];");
    AppendTo(outr, "\n];");
    CloseStream( outl );
    CloseStream( outr );
end;

CaratCatReadEncodedLeft := function()
    if not ( IsRecord(CDB) and IsBound(CDB.root_dir) ) then
        Error("CDB record not initialized. Run CaratCatInit() first.");
    fi;
    Info(InfoCaratCat, 1, "Reading (left) ZClasses (Bravais) data ...");
    Read( Concatenation( CDB.root_dir, "/z.left.grp" ) );
    Info(InfoCaratCat, 1, "Reading (left) Cocycle data ...");
    Read( Concatenation( CDB.root_dir, "/c.left.grp" ) );
    Info(InfoCaratCat, 1, "Reading Bieberbach group info ...");
    Read( Concatenation( CDB.root_dir, "/torsionfree.grp" ) );
end;

CaratCatReadEncodedRight := function()
    if not ( IsRecord(CDB) and IsBound(CDB.root_dir) ) then
        Error("CDB record not initialized. Run CaratCatInit() first.");
    fi;
    Info(InfoCaratCat, 1, "Reading (right) ZClasses (Bravais) data ...");
    Read( Concatenation( CDB.root_dir, "/z.right.grp" ) );
    Info(InfoCaratCat, 1, "Reading (right) Cocycle data ...");
    Read( Concatenation( CDB.root_dir, "/c.right.grp" ) );
    Info(InfoCaratCat, 1, "Reading Bieberbach group info ...");
    Read( Concatenation( CDB.root_dir, "/torsionfree.grp" ) );
end;